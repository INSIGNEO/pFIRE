#!/usr/bin/env python3

import os
import subprocess as sp

import numpy as np

from configobj import ConfigObj

import flannel.io as fio

from .image_routines import load_image


class ResultObject:
    """
    Small object to hold registration result info
    """

    def __init__(self, registered_path, map_path, log_path, fixed_path=None,
                 moved_path=None):
        self.registered_path = registered_path
        self.map_path = map_path
        self.log_path = log_path
        self.fixed_path = fixed_path
        self.moved_path = moved_path


class pFIRERunnerMixin:
    """ Mixin class to provide a pFIRE runner interface
    """

    def __init__(self, *args, **kwargs):
        super(pFIRERunnerMixin, self).__init__(*args, **kwargs)
        self.pfire_fixed_path = None
        self.pfire_moved_path = None
        self.pfire_mask_path = None
        self.pfire_reg_path = None
        self.pfire_map_path = None


    def run_pfire(self, config_path, comm_size=1):
        """ Run pFIRE using provided config file
        """
        if comm_size != 1:
            raise RuntimeError("MPI pFIRE runs not yet supported")

        pfire_workdir, pfire_config = [os.path.normpath(x) for x in
                                       os.path.split(config_path)]

        # opening explicitly causes failure on file nonexistence
        with open(config_path, 'r') as fh:
            config = ConfigObj(fh)
        print("Running pFIRE on {}".format(pfire_config))
      
        self.pfire_fixed_path = os.path.join(pfire_workdir, config['fixed'])
        self.pfire_moved_path = os.path.join(pfire_workdir, config['moved'])
             
        # TODO Check Registered image is in xdmf format
        # TODO Check registration Map is in xdmf format
        
        # Log file name
        self.pfire_logfile = os.path.join(
            pfire_workdir,
            "{}_pfire.log".format(os.path.splitext(pfire_config)[0]))
        
        # Run pFIRE
        with open(self.pfire_logfile, 'w') as logfile:
            pfire_args = ['pfire', pfire_config]
            res = sp.run(pfire_args, cwd=pfire_workdir, stdout=logfile,
                         stderr=logfile)
        if res.returncode != 0:
            raise RuntimeError("Failed to run pFIRE, check log for details: {}"
                               "".format(self.pfire_logfile))

        # Get pFIRE registered image path filename
        try:
            self.pfire_reg_path = os.path.join(pfire_workdir, config['registered'])
            print("self.pfire_reg_path="+ self.pfire_reg_path)
        except KeyError:
            pass
            
        # Get pFIRE registered map  path filename           
        try:
            self.pfire_map_path = os.path.join(pfire_workdir, config['map'])
            print("self.pfire_map_path="+ self.pfire_map_path)
        except KeyError:
            pass
           
        try:
            self.pfire_mask_path = os.path.join(pfire_workdir, config['mask'])
        except KeyError:
            pass

        if not (self.pfire_reg_path or self.pfire_map_path):
            raise RuntimeError("Failed to find registered image or map named in config file")
        

class ShIRTRunnerMixin:
    """ Mixin class to provide a ShIRT runner interface accepted a pFIRE config
    file
    """

    default_mask_name = "default_mask.mask"
    config_defaults = {"mask": None}

    def __init__(self, *args, **kwargs):
        super(ShIRTRunnerMixin, self).__init__(*args, **kwargs)
        self.shirt_fixed_path = None
        self.shirt_moved_path = None
        self.shirt_mask_path = None
        self.shirt_reg_path = None
        self.shirt_map_path = None


    def _build_shirt_config(self, config_file):
        """
        Build ShIRT config object given pfire config file
        """
        # read pfire config and merge with default options to get full option list
        defaults = ConfigObj(self.config_defaults)
        config = ConfigObj(config_file)
        defaults.merge(config)
        return defaults

    @staticmethod
    def _strip_imgname(name):
        shirtified = os.path.basename(name)
        return os.path.splitext(shirtified)[0]

    def run_shirt(self, config_path):
        """
        Run ShIRT using a pfire config file for input
        """
        config = self._build_shirt_config(config_path)
        work_dir, config_file = os.path.split(config_path)
        
        # FIXME workdire is empty if config path is relatove. make it absolute
        
        print("Running ShIRT on {}".format(config_file))

        self.shirt_reg_path = 'shirt_{}_registered.image'.format(self.name)
        self.shirt_map_path = 'shirt_{}_map.map'.format(self.name)

        for fom in ['fixed', 'moved', 'mask']:
            if fom == 'mask' and config[fom] is None:
                continue
            if config[fom].endswith(".image"):
                setattr(self, "shirt_{}_path".format(fom),
                        os.path.join(work_dir, config[fom]))
            else:
                data = load_image(os.path.join(work_dir, config[fom]))
                newname = os.path.join(
                    work_dir, os.path.splitext(config[fom])[0] + '.image')
                setattr(self, "shirt_{}_path".format(fom), newname)
                fio.save_image(data, newname)

        if config['mask'] is None:
            self.shirt_mask_path = os.path.join(work_dir,
                                                self.default_mask_name)
            data = load_image(os.path.join(work_dir, config['fixed']))
            data = np.full(data.shape, 1.0)
            fio.save_image(data, self.shirt_mask_path)

        shirt_args = ['ShIRT', 'Register', 'verbose',
                      'Fixed', self._strip_imgname(self.shirt_fixed_path),
                      'Moved', self._strip_imgname(self.shirt_moved_path),
                      'Mask', self._strip_imgname(self.shirt_mask_path),
                      'NodeSpacing', config['nodespacing'],
                      'Registered', self._strip_imgname(self.shirt_reg_path),
                      'Map', self._strip_imgname(self.shirt_map_path)]



        # ADd test FIXME
        #SHIRT_BIN
        
        shirt_env = {}
        SHIRT_DIR="/home/tartarini/LOCAL/app/Shirt/"#~/LOCAL/app/Shirt/bin/ShIRT

        #shirt_env['DATAPATH'] = '/home/tartarini/fixing_bugs_pfire/pFIRE/benchmarking/brain2d/'
        shirt_env['DATAPATH'] = './'
        shirt_env['DISPLAYPATH'] = shirt_env['DATAPATH'] +"display/"

        #SHIRT_DIR=/home/tartarini/ShIRT/
        shirt_env['SYSPATH']= SHIRT_DIR+"bin/"
        shirt_env['SCRIPTPATH']=SHIRT_DIR+"Scripts/"
        shirt_env['PATH']= shirt_env['SYSPATH']
        #export PATH=$PATH:$SHIRT_DIR/1.1/bin/


        self.shirt_logfile = os.path.join(
            work_dir, "{}_shirt.log".format(os.path.splitext(config_file)[0]))
        with open(self.shirt_logfile, 'w') as logfile:
            print("work_dir=" + work_dir)
            print(shirt_env)
            
            sp.run(['ShIRT', 'setpath', 'DataPath', '.'])#, cwd=work_dir)

            res = sp.run(shirt_args, stdout=logfile, stderr=logfile,  env=shirt_env)
                         #cwd=work_dir)

        if res.returncode != 0:
            raise RuntimeError("Failed to run ShIRT, check log for details: {}"
                               "".format(self.shirt_logfile))

        self.shirt_fixed_path = os.path.join(work_dir, self.shirt_fixed_path)
        self.shirt_moved_path = os.path.join(work_dir, self.shirt_moved_path)
        self.shirt_mask_path = os.path.join(work_dir, self.shirt_mask_path)
        self.shirt_reg_path = os.path.join(work_dir, self.shirt_reg_path)
        self.shirt_map_path = os.path.join(work_dir, self.shirt_map_path)
