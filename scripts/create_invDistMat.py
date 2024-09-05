import argparse
import pyreadr
import pandas as pd
import numpy as np
from scipy.spatial import distance_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import save_npz
from scipy.spatial import cKDTree

def get_arguments():
	parser = argparse.ArgumentParser(description="Create inverse distance matrix from cell spatial locations")
	parser.add_argument("-i", "--input", help="cell spatial location matrix in rds format", required=True, type=str)
	parser.add_argument("-s", "--scale", help="value to scale matrix by (5 recommended)", required=True, type=int)
	parser.add_argument("-d", "--distance", help="quantile of values to keep (0.05 recommended)", required=True, type=float)
	parser.add_argument("-o", "--outname", help="name of output file (in npz format)", required=True, type=str)
	return parser.parse_args()

args=get_arguments()

distance_cut = args.distance
scale = args.scale

result=pyreadr.read_r(args.input)
spatLocs=result[None]
spatLocs=spatLocs.set_index("cell_ID")
distMat=pd.DataFrame(distance_matrix(spatLocs.values, spatLocs.values), index=spatLocs.index, columns=spatLocs.index)
invDistMat=distMat.max().max() - distMat
cutoff=np.quantile(invDistMat, (1-distance_cut))
invDistMat=invDistMat.where(invDistMat>=cutoff, other=0)
invDistMat=invDistMat.div(scale)
np.fill_diagonal(invDistMat.values, 0)

sparseMat = csr_matrix(invDistMat)

save_npz(args.outname, sparseMat)
