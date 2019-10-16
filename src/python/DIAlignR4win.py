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
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
import math

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
			if file.endswith(".chrom.mzML"):
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

def extractXIC_group(mz, chromIndices):
	XIC_group = [None]*len(chromIndices)
	for i in range(len(chromIndices)):
		XIC_group[i] = mz.getChromatogram(chromIndices[i]).get_peaks() # Chromatogram is a tuple (time_array, intensity_array) of numpy array
	return XIC_group

def main():
	parser = ArgumentParser(description='Align XICs of precursors across multiple Targeted-MS runs and outputs quantitative data matrix.', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-in1', '--in_dataPath', metavar = '', help="location of parent directory of osw and mzml files.\n",  type = Path, default = None, required = True)
	parser.add_argument('-in2', '--alignType' , metavar = 'hybrid', type = str, nargs = 1, default = 'hybrid', required = False, choices=['hybrid', 'local', 'global'], help = 'Type of retention time alignment.\n')
	parser.add_argument('--maxFdrQuery' , metavar = '0.05', type = float, nargs = 1, default = 0.05, required = False, help = 'Type of retention time alignment.\n')
	parser.add_argument('--maxFdrLoess' , metavar = '0.01', type = float, nargs = 1, default = 0.01, required = False, help = 'Type of retention time alignment.\n')
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
	maxFdrLoess = args.maxFdrLoess
	spanvalue = args.spanvalue
	normalization = args.normalization
	simMeasure = args.simMeasure
	SgolayFiltOrd = args.SgolayFiltOrd
	SgolayFiltLen = args.SgolayFiltLen
	goFactor = args.goFactor
	geFactor = args.geFactor
	cosAngleThresh = args.cosAngleThresh
	OverlapAlignment = args.OverlapAlignment
	dotProdThresh = args.dotProdThresh
	gapQuantile = args.gapQuantile
	hardConstrain = args.hardConstrain
	samples4gradient = args.samples4gradient
	expRSE = args.expRSE
	samplingTime = args.samplingTime
	RSEdistFactor = args.RSEdistFactor

	oswFiles = []
	mzmlFiles = []
	runs = []
	peptides = []
	if args.in_dataPath.exists():
		curPath = args.in_dataPath
		runs = getRunNames(curPath, oswFiles, mzmlFiles)

	runs = ['170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#2_400-650mz_msms35', '170413_AM_BD-ZH12_BoneMarrow_W_10%_DIA_#1_400-650mz_msms34']
	# Get indices of chromatograms for each peptide.
	oswFiles = dict.fromkeys(runs)
	print("Reading osw and mzML files to fetch nativeIDs")
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
		mzml_db = os.path.join(curPath, 'mzml', run + '.chrom.mzML')
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
					chromatogramIndex[index] = str(i)
		df['chromatogramIndex'] = chromatogramIndex 
		df = df.assign(chromatogramIndex = lambda x: df.groupby(['transition_group_id', 'peak_group_rank'], sort=False)\
          .transform(lambda idx: ','.join(idx))['chromatogramIndex']).drop('transition_id', axis=1)
		transition_ids = df[df.m_score < 0.01]['transition_group_id'].tolist()
		peptides =  list(set(peptides) | set(df[df.m_score < 0.01]['transition_group_id'].tolist()))
		df.drop_duplicates(inplace=True)
		oswFiles[run] = df.reset_index(drop= True)
		print("Fetched data from {0}".format(run))
	peptides.sort()
	
	print("Collecting pointers to mzML files.")
	lowess_sm = sm.nonparametric.lowess
	mzPntrs = dict.fromkeys(runs)
	filter = SavitzkyGolayFilter()
	p = filter.getParameters()
	p.update({b'frame_length': args.SgolayFiltLen, b'polynomial_order': args.SgolayFiltOrd})
	filter.setParameters(p)
	for run in mzPntrs.keys():
		mzml_db = os.path.join(curPath, 'mzml', run + '.chrom.mzML')
		mz = OnDiscMSExperiment()
		mz.openFile(mzml_db)
		#filter.filterExperiment(mz)
		mzPntrs[run] = mz
	print("Pointers collected.")

	num_of_pep = len(peptides)
	rtTbl = {run: [None]*num_of_pep for run in runs}
	intesityTbl = {run: [None]*num_of_pep for run in runs}
	lwTbl = {run: [None]*num_of_pep for run in runs}
	rwTbl = {run: [None]*num_of_pep for run in runs}

	print("Aligning analytes using reference-based approach.")
	loessFits = dict()
	for pepIdx in range(2):
		peptide = peptides[pepIdx]
		# Select reference run based on m-score
		minMscore = 1.0
		minrunIdx = None
		for runIdx in range(len(runs)):
			df = oswFiles[runs[runIdx]]
			m_score = df.loc[(df.transition_group_id == peptide) & (df.peak_group_rank == 1), "m_score"].reset_index(drop= True)
			if not m_score.empty:
				m_score = m_score[0]
				if m_score < minMscore:
					minMscore = m_score
					minrunIdx = runIdx
		print(runs[minrunIdx])

		# Set the feature value as in the reference run
		ref = runs[minrunIdx]
		df_ref = oswFiles[ref]
		vec = df_ref.loc[(df_ref.transition_group_id == peptide) & (df_ref.peak_group_rank == 1), ['RT', 'Intensity', 'leftWidth', 'rightWidth']]
		vec = vec.reset_index(drop=True)
		rtTbl[ref][pepIdx] = vec['RT'].iloc[0]
		intesityTbl[ref][pepIdx] = vec['Intensity'].iloc[0]
		lwTbl[ref][pepIdx] = vec['leftWidth'].iloc[0]
		rwTbl[ref][pepIdx] = vec['rightWidth'].iloc[0]

		# Extract chromatograms from the reference run
		exps = set(runs) - set([ref])
		chromIndices = df_ref.loc[(df_ref.transition_group_id == peptide) & (df_ref.peak_group_rank == 1), 'chromatogramIndex'].to_string(header = False, index = False)
		chromIndices = [int(x) for x in chromIndices.split(",")]
		XICs_ref = extractXIC_group(mzPntrs[ref], chromIndices)
		# Align all runs to reference run
		for eXp in exps:
			# Get XIC_group from experiment run
			df_eXp = oswFiles[eXp]
			chromIndices = df_eXp.loc[(df_eXp.transition_group_id == peptide) & (df_eXp.peak_group_rank == 1), 'chromatogramIndex']
			if not chromIndices.empty:
				chromIndices = chromIndices.to_string(header = False, index = False)
				chromIndices = [int(x) for x in chromIndices.split(",")]
				XICs_eXp = extractXIC_group(mzPntrs[eXp], chromIndices)
				# Get the loess fit for hybrid alignment
				pair = '{0}_{1}'.format(ref, eXp) # Update these names with run0 run1 instead of filename. May be I can use series.
				if pair in loessFits.keys():
					Loess_fit = loessFits[pair]
				else:
					# Can't update dataframe after .loc in python?
					df_ref_rt = df_ref.loc[(df_ref.m_score < maxFdrLoess) & (df_ref.peak_group_rank == 1), ['transition_group_id', 'RT']]
					df_eXp_rt = df_eXp.loc[(df_eXp.m_score < maxFdrLoess) & (df_eXp.peak_group_rank == 1), ['transition_group_id', 'RT']]
					RUNS_RT = pd.merge(df_ref_rt, df_eXp_rt, how = 'inner', on = 'transition_group_id', suffixes=('_ref', '_eXp'))
					X = RUNS_RT.loc[:, 'RT_ref'].values.reshape(-1, 1)
					y = RUNS_RT.loc[:, 'RT_eXp'].values
					Loess_fit = LinearRegression().fit(X, y)
					loessFits[pair] = Loess_fit
				rse = min(9.0, expRSE)
				noBeef = math.ceil(RSEdistFactor*rse/samplingTime)
				tVec_ref = XICs_ref[0][0] # Extracting time component
				tVec_eXp = XICs_eXp[0][0] # Extracting time component
				B1p = Loess_fit.predict(np.array([[tVec_ref[0]]], dtype = 'float64'))
				B2p = Loess_fit.predict(np.array([[tVec_ref[-1]]], dtype = 'float64'))
				intensityList_ref = [XICs_ref[fragIdx][1] for fragIdx in range(len(XICs_ref))]
				intensityList_ref = np.ascontiguousarray(np.asarray(intensityList_ref, dtype = 'float64'))
				intensityList_eXp = [XICs_eXp[fragIdx][1] for fragIdx in range(len(XICs_eXp))]
				intensityList_eXp = np.ascontiguousarray(np.asarray(intensityList_eXp, dtype = 'float64'))
				obj = PyDIAlign.AffineAlignObj(500, 500)
				PyDIAlign.alignChromatogramsCppSimple(obj, intensityList_ref, intensityList_eXp,
	                                                  b'hybrid', tVec_ref, tVec_eXp, b'mean',
	                                                  b'dotProductMasked')
				AlignedIndices = pd.DataFrame(list(zip(obj.indexA_aligned, obj.indexB_aligned, obj.score)),
	                                          columns =['indexAligned_ref', 'indexAligned_eXp', 'score'])
							# C++ output is based on 1-index.
				tAligned_ref = mapIdxToTime(tVec_ref, AlignedIndices['indexAligned_ref']-1)
				tAligned_eXp = mapIdxToTime(tVec_eXp, AlignedIndices['indexAligned_eXp']-1)
				rtTbl[eXp][pepIdx] = tAligned_eXp[np.argmin(abs(tAligned_ref - rtTbl[ref][pepIdx]))]
				df = df_eXp.loc[(df_eXp.transition_group_id == peptide), ['RT', 'Intensity', 'leftWidth', 'rightWidth', 'peak_group_rank', 'm_score']]
				adaptiveRT = RSEdistFactor*rse
				df = df.loc[((abs(df.RT - rtTbl[eXp][pepIdx]) <= adaptiveRT) & (df.m_score < maxFdrQuery) & df.peak_group_rank.idxmin()), ]
				if not df.empty:
					lwTbl[eXp][pepIdx] = df.loc[0, 'leftWidth']
					rtTbl[eXp][pepIdx] = df.loc[0, 'RT']
					rwTbl[eXp][pepIdx] = df.loc[0, 'rightWidth']
					intesityTbl[eXp][pepIdx] = df.loc[0, 'Intensity']
	# Saving tables
	df = pd.DataFrame.from_dict(rtTbl, orient='columns')
	df['peptides'] = peptides
	df.set_index('peptides', inplace=True)
	df.to_csv('rtTblPy.csv', na_rep = 'NA')

	df = pd.DataFrame.from_dict(intesityTbl, orient='columns')
	df['peptides'] = peptides
	df.set_index('peptides', inplace=True)
	df.to_csv('intesityTblPy.csv', na_rep = 'NA')

	df = pd.DataFrame.from_dict(lwTbl, orient='columns')
	df['peptides'] = peptides
	df.set_index('peptides', inplace=True)
	df.to_csv('lwTblPy.csv', na_rep = 'NA')

	df = pd.DataFrame.from_dict(rwTbl, orient='columns')
	df['peptides'] = peptides
	df.set_index('peptides', inplace=True)
	df.to_csv('rwTblPy.csv', na_rep = 'NA')
	return 1


if __name__=='__main__':
	main()
