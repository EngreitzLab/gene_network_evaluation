import os
import argparse

# Rename function to compute_{eval_measure}
def main(mdata, optional=None, 
	 prog_key='prog', rna_key='rna', atac_key='atac', 
	 inplace=True):
    
    #TODO: Don't copy entire mudata only relevant Dataframe
    mdata = mdata.copy() if not inplace else mdata
  
    # Compute & update mudata with eval measure
    mdata[prog_key].var['eval_measure'] = None

    if not inplace: return mdata[prog_key].var['eval_measure']
	
if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('mudataObj_path')
    parser.add_argument('-o', '--optional') 
    parser.add_argument('-pk', '--prog_key', default='prog', type=str) 
    parser.add_argument('-rk', '--rna_key', default='rna', type=str) 
    parser.add_argument('-ak', '--atac_key', default='atac', type=str) 
    parser.add_argument('--output', action='store_false') 

    args = parser.parse_args()

    import mudata
    mdata = mudata.read(args.mudataObj_path)

    main(mdata, optional=args.optional,
         prog_key=args.prog_key, rna_key=args.rna_key, 
         atac_key=args.atac_key, inplace=args.output)

