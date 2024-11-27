from Xdmf import *

if __name__ == "__main__":

    fileName ="TestXdmfRead.xmf"

    # create a simple empty file
    domain = XdmfDomain.New()
    collection = XdmfGridCollection.New()
    grid = XdmfUnstructuredGrid.New()
    attribute1 = XdmfAttribute.New()
    attribute2 = XdmfAttribute.New()
    information = XdmfInformation.New()

    domain.insert(collection)
    collection.insert(grid)
    grid.insert(attribute1)
    grid.insert(attribute2)
    grid.insert(information)
    
    writer = XdmfWriter.New(fileName)
    domain.accept(writer)

    # read file using XPaths and verify downcasts to appropriate XdmfItems
    reader = XdmfReader.New()

    domain = reader.read(fileName, "/Xdmf/Domain")
    print(str(len(domain)) + " ?= " + str(1))
    print("? " + str(isinstance(domain[0], XdmfDomain)))
    assert(len(domain) == 1)
    assert(isinstance(domain[0], XdmfDomain))

    collection = reader.read(fileName, "/Xdmf/Domain/Grid")
    print(str(len(collection)) + " ?= " + str(1))
    print("? " + str(isinstance(collection[0], XdmfGridCollection)))
    assert(len(collection) == 1)
    assert(isinstance(collection[0], XdmfGridCollection))
    
    grid = reader.read(fileName, "/Xdmf/Domain/Grid/Grid")
    print(str(len(grid)) + " ?= " + str(1))
    print("? " + str(isinstance(grid[0], XdmfUnstructuredGrid)))
    assert(len(grid) == 1)
    assert(isinstance(grid[0], XdmfUnstructuredGrid))

    attributes = reader.read(fileName, "/Xdmf/Domain/Grid/Grid/Attribute")
    print(str(len(attributes)) + " ?= " + str(2))
    print("? " + str(isinstance(attributes[0], XdmfAttribute)))
    print("? " + str(isinstance(attributes[1], XdmfAttribute)))
    assert(len(attributes) == 2)
    assert(isinstance(attributes[0], XdmfAttribute))
    assert(isinstance(attributes[1], XdmfAttribute))

    information = reader.read(fileName, "/Xdmf/Domain/Grid/Grid/Information")
    print(str(len(information)) + " ?= " + str(1))
    print("? " + str(isinstance(information[0], XdmfInformation)))
    assert(len(information) == 1)
    assert(isinstance(information[0], XdmfInformation))

