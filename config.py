"""config.py - Configuration module for SYGNAL
"""
import argparse
import json
import os
import shutil

class SygnalConfig:
    def __init__(self, config):
        self.config = config
        self.basedir = os.path.dirname(os.getcwd()) + '/'
        self.outdir = os.path.abspath(config['outdir'])
        self.tmpdir = os.path.abspath(config['tmpdir'])

    def outdir_path(self, path):
        return os.path.join(self.outdir, path)

    def tmpdir_path(self, path):
        return os.path.join(self.tmpdir, path)

    def clear_tmp(self):
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)

    def create_tmpdir(self, path):
        tmppath = os.path.join(self.tmpdir, path)
        if not os.path.exists(tmppath):
            os.makedirs(tmppath)



def sygnal_init():
    """Call this at the start of the program.
    Returns a SygnalConfig object"""
    parser = argparse.ArgumentParser(description="sygnal.py - Systems Genetic Network AnaLyse pipeline")
    parser.add_argument('--config',  default="sygnal_config.json", help="config file")
    parser.add_argument('--out',  default='output', help='output directory')
    parser.add_argument('--tmp',  default='tmp', help='temporary output directory')
    args = parser.parse_args()

    with open(args.config, 'r') as infile:
        config = json.load(infile)
        config['outdir'] = args.out
        config['tmpdir'] = args.tmp

        if not os.path.exists(config['outdir']):
            os.makedirs(config['outdir'])
        return SygnalConfig(config)
    
