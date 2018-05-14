import os
import errno

class cd:
   """Context manager for changing the current working directory"""
   def __init__(self, newPath):
      self.make_sure_path_exists(newPath)
      self.newPath = os.path.expanduser(newPath)

   def __enter__(self):
      self.savedPath = os.getcwd()
      os.chdir(self.newPath)

   def __exit__(self, etype, value, traceback):
      os.chdir(self.savedPath)

   def make_sure_path_exists(self,path):
      try:
         os.makedirs(path)
      except OSError as exception:
         if exception.errno != errno.EEXIST:
            raise
