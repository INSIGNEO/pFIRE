//
//   Copyright 2019 University of Sheffield
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#include "xdmfwriter.hpp"

#include <algorithm>
#include <sstream>

#include <boost/filesystem.hpp>

#include "image.hpp"
#include "indexing.hpp"
#include "infix_iterator.hpp"
#include "map.hpp"

namespace bf = boost::filesystem;

const std::string XDMFWriter::writer_name = "xdmf";
const std::vector<std::string> XDMFWriter::extensions = {".xdmf"};

XDMFWriter::XDMFWriter(std::string filespec, const MPI_Comm &comm)
  : HDFWriter(h5name_from_xdmfname(filespec), comm), rank(-1),
    xdmf_filename(split_filespec(filespec).first), xdmf_groupname(split_filespec(filespec).second),
    xdmf_tree(pt::ptree())
{
  // Only want to do reading/writing of the text file from rank 0
  // but still to hdf5 writing from all ranks
  MPI_Comm_rank(comm, &rank);

  // Pointless if rank != 0
  if (rank == 0)
  {
    prepopulate_xdmf_ptree();
  }
}

XDMFWriter::~XDMFWriter()
{
}

std::string XDMFWriter::write_image(const Image &image)
{
  // First need to write the hdf data
  HDFWriter::write_image(image);

  if (xdmf_groupname == "")
  {
    xdmf_groupname = "image";
  }

  // Assuming success, if rank 0 go ahead and populate grid information
  if (rank == 0)
  {
    pt::ptree &domain_tree = xdmf_tree.get_child("Xdmf.Domain");
    pt::ptree &grid = domain_tree.add("Grid", "");
    grid.add("<xmlattr>.Name", xdmf_groupname + std::string("_grid"));
    grid.add("<xmlattr>.GridType", "Uniform");

    // Define topology type
    {
      pt::ptree &topo = grid.add("Topology", "");
      topo.add("<xmlattr>.TopologyType", "3DCoRectMesh");
      intvector meshdims(image.shape());
      std::for_each(meshdims.begin(), meshdims.end(), [](integer &i) { i += 1; });
      std::ostringstream topodims;
      std::copy(
          meshdims.cbegin(), meshdims.cend(), infix_ostream_iterator<integer>(topodims, " "));
      topo.add("<xmlattr>.Dimensions", topodims.str());
    }

    // Define geometry
    {
      pt::ptree &geom = grid.add("Geometry", "");
      geom.add("<xmlattr>.GeometryType", "Origin_DxDyDz");
      // Origin
      pt::ptree &origin = geom.add("DataItem", "0 0 0");
      origin.add("<xmlattr>.Name", "Origin");
      origin.add("<xmlattr>.Dimensions", "3");
      origin.add("<xmlattr>.NumberType", "Float");
      origin.add("<xmlattr>.Format", "XML");
      // Spacing
      pt::ptree &spacing = geom.add("DataItem", "1 1 1");
      spacing.add("<xmlattr>.Name", "Spacing");
      spacing.add("<xmlattr>.Dimensions", "3");
      spacing.add("<xmlattr>.NumberType", "Float");
      spacing.add("<xmlattr>.Format", "XML");
    }

    // Finally add image data attribute
    {
      pt::ptree &attr = grid.add("Attribute", "");
      attr.add("<xmlattr>.Name", xdmf_groupname);
      attr.add("<xmlattr>.AttributeType", "Scalar");
      attr.add("<xmlattr>.Center", "Cell");

      // And link to hdf5 data
      std::ostringstream datapath;
      datapath << bf::path(h5_filename).filename().string() << ":" << h5_groupname;
      std::ostringstream imgdims;
      std::copy(image.shape().cbegin(), image.shape().cend(),
          infix_ostream_iterator<integer>(imgdims, " "));

      pt::ptree &imdata = attr.add("DataItem", datapath.str());
      imdata.add("<xmlattr>.Format", "HDF");
      imdata.add("<xmlattr>.Dimensions", imgdims.str());
    }

    // potentially fragile undocumented settings
    pt::xml_writer_settings<std::string> settings(' ', 4);
    pt::write_xml(xdmf_filename, xdmf_tree, std::locale(), settings);
  }

  std::string filepath = xdmf_filename + ":" + h5_groupname;
  return filepath;
}

std::string XDMFWriter::write_map(const Map &map)
{
  HDFWriter::write_map(map);

  if (xdmf_groupname == "")
  {
    xdmf_groupname = "map";
  }

  // Assuming success, if rank 0 go ahead and populate grid information
  if (rank == 0)
  {
    pt::ptree &domain_tree = xdmf_tree.get_child("Xdmf.Domain");
    pt::ptree &grid = domain_tree.add("Grid", "");
    grid.add("<xmlattr>.Name", xdmf_groupname + std::string("_grid"));
    grid.add("<xmlattr>.GridType", "Uniform");

    // Define topology type
    {
      pt::ptree &topo = grid.add("Topology", "");
      topo.add("<xmlattr>.TopologyType", "3DCoRectMesh");
      intvector meshdims(map.shape());
      std::ostringstream topodims;
      std::copy(
          meshdims.cbegin(), meshdims.cend(), infix_ostream_iterator<integer>(topodims, " "));
      topo.add("<xmlattr>.Dimensions", topodims.str());
    }

    // Define geometry
    {
      pt::ptree &geom = grid.add("Geometry", "");
      geom.add("<xmlattr>.GeometryType", "Origin_DxDyDz");
      // Origin
      std::ostringstream originss;
      floatvector originvec = map.low_corner();
      std::copy(
          originvec.cbegin(), originvec.cend(), infix_ostream_iterator<floating>(originss, " "));
      pt::ptree &origin = geom.add("DataItem", originss.str());
      origin.add("<xmlattr>.Name", "Origin");
      origin.add("<xmlattr>.Dimensions", "3");
      origin.add("<xmlattr>.NumberType", "Float");
      origin.add("<xmlattr>.Format", "XML");
      // Spacing
      std::ostringstream spacingss;
      std::copy(map.spacing().cbegin(), map.spacing().cend(),
          infix_ostream_iterator<floating>(spacingss, " "));
      pt::ptree &spacing = geom.add("DataItem", spacingss.str());
      spacing.add("<xmlattr>.Name", "Spacing");
      spacing.add("<xmlattr>.Dimensions", "3");
      spacing.add("<xmlattr>.NumberType", "Float");
      spacing.add("<xmlattr>.Format", "XML");
    }

    // Finally add map data attribute
    {
      pt::ptree &attr = grid.add("Attribute", "");
      attr.add("<xmlattr>.Name", xdmf_groupname);
      attr.add("<xmlattr>.AttributeType", "Vector");
      attr.add("<xmlattr>.Center", "Node");

      // Have separate datasets in hdf5 so use join
      std::ostringstream joinfunc;
      joinfunc << "join(";
      for (uinteger idx = 0; idx < map.ndim(); idx++)
      {
        if (idx != 0)
        {
          joinfunc << ", ";
        }
        joinfunc << "$"
                 << static_cast<unsigned long>(
                        idx); // resolve ambiguity due to lack of explicit << overload
      }
      joinfunc << ")";
      std::ostringstream mapdims;
      std::copy(
          map.shape().cbegin(), map.shape().cend(), infix_ostream_iterator<integer>(mapdims, " "));
      std::ostringstream vecdims;
      vecdims << mapdims.str() << " " << static_cast<unsigned>(map.ndim());
      pt::ptree &vecdat = attr.add("DataItem", "");
      vecdat.add("<xmlattr>.ItemType", "Function");
      vecdat.add("<xmlattr>.Dimensions", vecdims.str());
      vecdat.add("<xmlattr>.Function", joinfunc.str());

      // And link to hdf5 datasets
      for (uinteger idx = 0; idx < map.ndim(); idx++)
      {
        std::ostringstream datapath;
        datapath << bf::path(h5_filename).filename().string() << ":" << h5_groupname << "/" << _components[idx];

        pt::ptree &imdata = attr.add("DataItem", datapath.str());
        imdata.add("<xmlattr>.Name", std::string("d") + _components[idx]);
        imdata.add("<xmlattr>.Format", "HDF");
        imdata.add("<xmlattr>.Dimensions", mapdims.str());
      }
    }

    // potentially fragile undocumented settings
    pt::xml_writer_settings<std::string> settings(' ', 4);
    pt::write_xml(xdmf_filename, xdmf_tree, std::locale(), settings);
  }

  std::string filepath = xdmf_filename + ":" + h5_groupname;
  return filepath;
}

std::string XDMFWriter::h5name_from_xdmfname(const std::string &filename)
{
  string_pair pp = split_filespec(filename);
  return pp.first + std::string(".h5") + ":" + pp.second;
}

void XDMFWriter::prepopulate_xdmf_ptree()
{
  pt::ptree domain_tree;
  xdmf_tree.put_child("Xdmf.Domain", domain_tree);
  xdmf_tree.add("Xdmf.<xmlattr>.Version", "3.0");
}
