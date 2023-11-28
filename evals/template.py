import os
import argparse

def main(adata, optional=None, inplace=True):
  
  adata = adata.copy() if inplace else adata
  
  # Update adata with eval measure
  adata.var['eval_measure'] = None

  if inplace: return adata.var['eval_measure']
	
if __name__=='__main__':
  parser = argparse.ArgumentParser()

	parser.add_argument('anndataObj')
	parser.add_argument('-o', '--optional') 

	args = parser.parse_args()
    
	main(args.anndataObj, optional=args.optional)
