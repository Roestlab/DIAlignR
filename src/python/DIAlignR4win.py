#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 13:14:16 2018

@author: shubham
"""
import os
import sqlite3
from sqlite3 import Error as sql_error
import pandas as pd
from pyopenms import *
import scipy
from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path
import datetime

def getRunNames(dataPath, oswFiles, mzmlFiles):
	runs = []
	oswdataPath = os.path.join(dataPath, 'osw')
	mzmldataPath = os.path.join(dataPath, 'mzml')
	for root, dirs, files in os.walk(oswdataPath):
		for file in files:
			if file.endswith(".osw"):
				oswFiles.append(os.path.join(root, file))
				runs.append(file[:-4])

	runs2 = []
	for root, dirs, files in os.walk(mzmldataPath):
		for file in files:
			if file.endswith("_chrom.mzML"):
				mzmlFiles.append(os.path.join(root, file))
				runs2.append(file[:-11])
	if set(runs) == set(runs2):
		print("All files are matched")
	return(runs)

def getQuery(maxFdrQuery, peptides = None):
	if peptides is None:
		selectPeptide  = ''
	else:
		selectPeptide = " AND transition_group_id IN ({})".format(peptides)

	query = """SELECT PEPTIDE.MODIFIED_SEQUENCE || '_' || PRECURSOR.CHARGE AS transition_group_id,
	RUN.FILENAME AS filename,
	FEATURE.EXP_RT AS RT,
	FEATURE.DELTA_RT AS delta_rt,
	PRECURSOR.LIBRARY_RT AS assay_RT,
	FEATURE_MS2.AREA_INTENSITY AS Intensity,
	FEATURE.LEFT_WIDTH AS leftWidth,
	FEATURE.RIGHT_WIDTH AS rightWidth,
	SCORE_MS2.RANK AS peak_group_rank,
	SCORE_MS2.QVALUE AS m_score,
	TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AS transition_id
	FROM PRECURSOR
	INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID AND PRECURSOR.DECOY=0
	INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID
	INNER JOIN FEATURE ON FEATURE.PRECURSOR_ID = PRECURSOR.ID
	INNER JOIN RUN ON RUN.ID = FEATURE.RUN_ID
	INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID = PRECURSOR.ID
	LEFT JOIN FEATURE_MS2 ON FEATURE_MS2.FEATURE_ID = FEATURE.ID
	LEFT JOIN SCORE_MS2 ON SCORE_MS2.FEATURE_ID = FEATURE.ID
	WHERE SCORE_MS2.QVALUE < {0:f} {1:s}
	ORDER BY transition_group_id,
	peak_group_rank;""".format(maxFdrQuery, selectPeptide)


	return(query)

def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except sql_error as e:
        print(e)
    return conn

def loadmzMLdata(filename, mz):
    print("Succeded in loading {0:s}".format(filename))
    return True

def SGolayOnRuns(runs):
    filter = SavitzkyGolayFilter()
    p = filter.getParameters()
    p.update({'frame_length': 9})
    filter.setParameters(p)
    # filter.getParameters().getValue('frame_length')
    for run in runs:
        filter.filterExperiment(run)
    return True
    

def SummaryOfRuns(runs):
    print("List of runs = " + str(runs))
    for i in range(len(runs)):
        print ("Number of Chromatograms in run[" + str(i) + "] = " + str(runs[i].getNrChromatograms())) # 140080, 140080
    print("")
    return True


def main():
	parser = ArgumentParser(description='Align XICs of precursors across multiple Targeted-MS runs and outputs quantitative data matrix.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-in1', '--in_dataPath', metavar = '', help="location of parent directory of osw and mzml files.\n",  type = Path, default = [None], required = True)
	parser.add_argument('-in2', '--alignType' , metavar = 'hybrid', type = str, nargs = 1, default = 'hybrid', required = False, choices=['hybrid', 'local', 'global'], help = 'Type of retention time alignment.\n')
	parser.add_argument('--maxFdrQuery' , metavar = '0.05', type = float, nargs = 1, default = 0.05, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--query' , metavar = 'query', type = str, nargs = 1, default = None, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--spanvalue' , metavar = '0.1', type = float, nargs = 1, default = 0.1, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--normalization' , metavar = 'mean', type = str, nargs = 1, default = 'mean', required = False, choices=['mean', 'l2', 'none'], help = 'Type of retention time alignment.\n')
	parser.add_argument('--simMeasure' , metavar = 'dotProductMasked', type = str, nargs = 1, default = 'dotProductMasked', required = False, choices=['dotProductMasked', 'dotProduct', 'cosineAngle',
																																						'cosine2Angle', 'euclideanDist', 'covariance', 'correlation'],
																																					 help = 'Method to build similarity matrix through XIC intensities.\n')
	parser.add_argument('--SgolayFiltOrd' , metavar = '4', type = int, nargs = 1, default = 4, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--SgolayFiltLen' , metavar = '9', type = int, nargs = 1, default = 9, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--goFactor' , metavar = '0.125', type = float, nargs = 1, default = 0.125, required = False, help = 'Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.\n')
	parser.add_argument('--geFactor' , metavar = '40', type = float, nargs = 1, default = 40, required = False, help = 'Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.\n')
	parser.add_argument('--cosAngleThresh' , metavar = '0.3', type = float, nargs = 1, default = 0.3, required = False, help = 'In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.\n')
	parser.add_argument('--OverlapAlignment' , metavar = 'True', type = bool, nargs = 1, default = True, required = False, help = 'An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.\n')
	parser.add_argument('--dotProdThresh' , metavar = '0.96', type = float, nargs = 1, default = 0.96, required = False, help = ' In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.\n')
	parser.add_argument('--gapQuantile' , metavar = '0.5', type = float, nargs = 1, default = 0.5, required = False, help = 'Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.\n')
	parser.add_argument('--hardConstrain' , metavar = 'False', type = bool, nargs = 1, default = False, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--samples4gradient' , metavar = '100', type = int, nargs = 1, default = 100, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--expRSE' , metavar = '8.0', type = float, nargs = 1, default = 8.0, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--samplingTime' , metavar = '3.4', type = float, nargs = 1, default = 3.4, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--RSEdistFactor' , metavar = '3.5', type = float, nargs = 1, default = 3.5, required = False, help = 'Type of retention time alignment.\n')
    
	args = parser.parse_args()
	maxFdrQuery = args.maxFdrQuery
	oswFiles = []
	mzmlFiles = []
	runs = []
	if args.in_dataPath.exists():
		curPath = args.in_dataPath
		runs = getRunNames(curPath, oswFiles, mzmlFiles)
	for i in range(len(runs)):
		print(runs[i])

	runs = ["170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35"]
	for run in runs:
		osw_db = os.path.join(curPath, 'osw', run + '.osw')
		# create a osw_db connection
		conn = create_connection(osw_db)
		cur = conn.cursor()
		query = getQuery(maxFdrQuery)
		cur.execute(query)
		rows = cur.fetchall() # Number of rows = 3338*5=16690
		conn.close() # close the osw_db connection
		# Convert list of tuples to dataframes.
		df = pd.DataFrame(rows, columns=['transition_group_id', 'filename', 'RT', 'delta_rt',
		 'assay_RT', 'Intensity', 'leftWidth', 'rightWidth', 'peak_group_rank', 'm_score', 'transition_id'])
		# Get ChromatogramID (native) and ChromatogramIndex
		mzml_db = os.path.join(curPath, 'mzml', run + '_chrom.mzML')
		mz = OnDiscMSExperiment()
		mz.openFile(mzml_db)
		meta_data = mz.getMetaData()
		# Find matching indices of chromatogram and add them to data-frame
		nrow = len(df.index)
		chromatogramIndex = [None]*nrow
		transition_ids = df['transition_id'].tolist()
		for i in range(meta_data.getNrChromatograms()):
			nativeID = int(meta_data.getChromatogram(i).getNativeID())
			if nativeID in transition_ids:
				indices = [index for index, value in enumerate(transition_ids) if value == nativeID]
				for index in indices:
					chromatogramIndex[index] = i
		df['chromatogramIndex'] = chromatogramIndex 
	
	for i in range(10):
		print(df[['transition_id', 'chromatogramIndex']].iloc[i])

if __name__=='__main__':
	main()



"""

filter = SavitzkyGolayFilter()
p = filter.getParameters()
# p.update({b'frame_length' : 9})
p.update({'frame_length' : 9})
filter.setParameters(p)
# filter.getParameters().getValue('frame_length')
for run in runs:
	filter.filterExperiment(run)
"""