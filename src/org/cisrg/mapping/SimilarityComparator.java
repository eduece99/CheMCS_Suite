package org.cisrg.mapping;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.cisrg.ambit.SmartsHelper;
import org.cisrg.hyperstructures.CDKSMARTSHyperstructureFitness;
import org.cisrg.knime.GraphSimilarityNodeModel.SimilarityType;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.smsd.AtomAtomMapping;

public class SimilarityComparator {
	
	

	/**
	 * 
	 */
	//private GraphSimilarityNodeModel graphSimilarityNodeModel;
	public SimilarityComparator( MCSMethods m, boolean weights, boolean ghostSubstructures, boolean detailedGhostInfo ) {
		
		//this.comparator = e;
		this.mcsMapper = null;
		this.useWeights = weights;
		this.ghostSubstructures = ghostSubstructures;
		
		this.uniqueSubstructureCount = 0;
		this.mcsGhostCount = 0;
		
		/*
		 * this.mappingAlgorithmName = e.getAlgorithmType(); this.useRaymondHeuristics =
		 * e.isCliqueDetectionRaymondHeuristics(); this.useRingHeuristics =
		 * e.isCliqueDetectionRingHeuristics(); this.topologicalDistanceLimit =
		 * e.getTopologicalDistanceLimit(); this.expansionTimeOut =
		 * e.getExpansionTimeOut();
		 */
		
		this.detailedGhostInfo = detailedGhostInfo;
		
		this.smaH = new SmartsHelper( DefaultChemObjectBuilder.getInstance() );
        this.sGenerator = new SmilesGenerator(SmiFlavor.Generic);
	}
	
	public SimilarityComparator( 
			MCSMethods m, 
			boolean weights, 
			boolean ghostSubstructures ,
			ExtendedAlgorithm mappingAlgorithmName,
			boolean useRaymondHeuristics,
			boolean useRingHeuristics,
			int topologicalDistanceLimit,
			int expansionTimeOut,
			boolean detailedGhostInfo
		) {
		
		this( m, weights, ghostSubstructures, mappingAlgorithmName, 
				useRaymondHeuristics, useRingHeuristics, topologicalDistanceLimit,
				expansionTimeOut, detailedGhostInfo, null, null
		);
		
		this.smaH = new SmartsHelper( DefaultChemObjectBuilder.getInstance() );
        this.sGenerator = new SmilesGenerator(SmiFlavor.Generic);
	}
	
	public SimilarityComparator( 
			MCSMethods m, 
			boolean weights, 
			boolean ghostSubstructures ,
			ExtendedAlgorithm mappingAlgorithmName,
			boolean useRaymondHeuristics,
			boolean useRingHeuristics,
			int topologicalDistanceLimit,
			int expansionTimeOut,
			boolean detailedGhostInfo,
			SmartsHelper sh,
			SmilesGenerator sg
		) {
		
	
		
		this.mcsMapper = m;
		this.useWeights = weights;
		this.ghostSubstructures = ghostSubstructures;
		
		this.uniqueSubstructureCount = 0;
		this.mcsGhostCount = 0;
		
		this.mappingAlgorithmName =  mappingAlgorithmName;
		this.useRaymondHeuristics = useRaymondHeuristics;
		this.useRingHeuristics = useRingHeuristics;
		this.topologicalDistanceLimit = topologicalDistanceLimit; 
		this.expansionTimeOut = expansionTimeOut;
		
		this.detailedGhostInfo = detailedGhostInfo;
		
		this.smaH = sh;
		this.sGenerator = sg;
	}
	
	
	
	public void calculateSimilarity( IAtomContainer rMol, IAtomContainer dbMol ) {
		
		if( mcsMapper != null ) {
			
			if( mcsMapper.getMainMol() != rMol )
				mcsMapper.setMainMol(rMol);
			
			mcsMapper.setQueryMol(dbMol);
			
			mcsMapper.execute();
			bondMaps = mcsMapper.getBestBondMatches();
			
			//List<List<Integer>> atomIndexLists = mcsMapper.getBestAtomIndexMatches();
			//atomMaps = new ArrayList<Map<IAtom, IAtom>>();  // FIXME
			for( List<Integer> indexList : mcsMapper.getBestAtomIndexMatches() ) {
				Map<IAtom, IAtom> atomMap = new HashMap<IAtom, IAtom>();
				
				for( int i = 0; i < indexList.size(); i++ ) {
					atomMap.put( rMol.getAtom(i) , dbMol.getAtom( indexList.get(i) ) );
				}

				//atomMaps.add(atomMap);  // FIXME
			}
			

			
			if( bondMaps == null || bondMaps.size() == 0 ) {
				tversky = 0.0;
				tanimoto = 0.0;
				
				return;
			}
			
			
			fragmentSizes = mcsMapper.fragmentSizes;
			mcsSMARTS =  mcsMapper.mcsSMARTS;
	
		} else {
			try {
				/*if( rMol instanceof IQueryAtomContainer ) {
					IQueryAtomContainer rMolQ = (IQueryAtomContainer) rMol;
					comparator.init(rMolQ, dbMol);
					logger.debug("IQueryAtomContainer used for similarity");
					//System.out.println("IQueryAtomContainer used for similarity");
				} else {
					comparator.init(rMol, dbMol, true, false);
					//System.out.println( "db atoms - " + dbMol.getAtomCount() );
				}*/
				
				
				ExtendedIsomorphism comparator = new ExtendedIsomorphism(
						rMol, 
						dbMol, 
						this.mappingAlgorithmName, 
						true, 
						false, 
						true, 
						this.useRaymondHeuristics,
						this.useRingHeuristics,
						this.topologicalDistanceLimit, 
						this.expansionTimeOut);
				//comparator.setChemFilters(true, true, true);
				
				bondMaps = comparator.getAllBondMaps();
				atomMaps = comparator.getAllAtomMapping();
				atomIndexMaps = comparator.getAllAtomIndexMapping();
				bondIndexMaps = comparator.getAllBondIndexMaps();
				
				fragmentSizes = comparator.getFragmentSizes();
				mcsSMARTS =  comparator.getMCSSMARTS();
				
				modProdConstructionTime = comparator.getModularProductConstructionTime();
				modProdNodeCount = comparator.getModularProductNodeCount();
				modProdEdgeDensity = comparator.getModularProductEdgeDensity();
				mcsExecTime = comparator.getElapsedTime();
			

			} catch (NullPointerException e3) {
				e3.printStackTrace();
				//logger.warn("Null pointer - similarity set to -1.0 " + rMol.getBondCount() + " " + dbMol.getBondCount() );
				//logger.warn( "common bonds - " + comparator.getAllBondMaps() );
				//logger.warn("settings - " + this.graphSimilarityNodeModel.mappingAlgorithmName  + " " + this.graphSimilarityNodeModel.useRaymondHeuristics + " " + this.graphSimilarityNodeModel.topologicalDistanceLimit + " " + this.graphSimilarityNodeModel.expansionTimeOut  );
				//logger.warn("stats - " + fragmentSizes  + " " + mcsSMARTS + " " + modProdConstructionTime + " " + modProdNodeCount + " " + modProdEdgeDensity + " " + mcsExecTime  );
				
				//similarity = 0.0;
			}
		}
		
		
		if( ! bondMaps.isEmpty() ) {
			
			try {
				
				IBond[] refBonds = bondMaps.get(0).keySet().toArray( new IBond[0] );
				IBond[] dbBonds = bondMaps.get(0).values().toArray( new IBond[0] );
				
				if( useWeights ) {
					
					commonBondWeights = new int[ refBonds.length ];
					commonBondTopologies = new String[ refBonds.length ];
					double weightBondSum = 0.0;
					double maxWeight = 0.0;  // we are not allowing negative weights
					
					for( int w = 0; w < bondMaps.get(0).size(); w++ ) {
						commonBondWeights[w] = (Integer) refBonds[w].getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType );
						commonBondTopologies[w] = "" + refBonds[w].getProperty( CDKSMARTSHyperstructureFitness.topologyType );
						weightBondSum += commonBondWeights[w];
						
						if( commonBondWeights[w] > maxWeight )
							maxWeight = commonBondWeights[w];
					}
					
					double weightRatio = ( weightBondSum / commonBondWeights.length ) / maxWeight;
					
					//tversky = ConvenienceTools.calculateWeightedTversky(weightBondSum, rMol.getBondCount(), dbMol.getBondCount(), 0.1, 0.9 );
					tversky = ConvenienceTools.calculateTversky( bondMaps.get(0).size(), rMol.getBondCount(), dbMol.getBondCount(), 0.1, 0.9 ) * weightRatio;
    				
					tanimoto = ConvenienceTools.calculateTversky( bondMaps.get(0).size(), rMol.getBondCount(), dbMol.getBondCount(), 1.0, 1.0 ) * weightRatio;
				} else {
					tversky = ConvenienceTools.calculateTversky( bondMaps.get(0).size(), rMol.getBondCount(), dbMol.getBondCount(), 0.1, 0.9 );
				
					tanimoto = ConvenienceTools.calculateTversky( bondMaps.get(0).size(), rMol.getBondCount(), dbMol.getBondCount(), 1.0, 1.0 );
					
					commonBondWeights = new int[ bondMaps.get(0).size() ];
					commonBondTopologies = new String[ bondMaps.get(0).size() ];
					for( int w = 0; w < bondMaps.get(0).size(); w++ ) {
						commonBondWeights[w] = 1;
						
						
						if( refBonds[w].getFlag( CDKConstants.ISINRING ) && dbBonds[w].getFlag( CDKConstants.ISINRING ) ) {
							commonBondTopologies[w] = "r";
						} else if( ! refBonds[w].getFlag( CDKConstants.ISINRING ) && ! dbBonds[w].getFlag( CDKConstants.ISINRING ) ) {
							commonBondTopologies[w] = "c";
						} else {
							commonBondTopologies[w] = "-";
						}
					}
				}
				
				
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		
		
    		if( ghostSubstructures ) {
    			
    			//IAtomContainer mcs = ConvenienceTools.createCommonSubgraph( dbMol, rMol, bondMaps.get(0) );
    			IAtomContainer mcs = ConvenienceTools.createCommonSubgraph( rMol, dbMol, bondMaps.get(0) );
    			
    			
    			if( this.detailedGhostInfo )
    				mcsGhostSMARTSSet = new ArrayList<String>( ghostRadius * mcs.getAtomCount() );
    			
    			//System.out.println( "hs mcs ghosts" );
    			
    			HashMap<Integer, Boolean> uniqueSubstructures = new HashMap<Integer, Boolean>();

    			for( int n=0; n < mcs.getAtomCount(); n++ ) {

    				IAtomContainer[] radii = ConvenienceTools.getNeighbourhoodGraphEdgeInduced(mcs, n, ghostRadius);
    				
    				for( int r = 1; r < radii.length; r++ ) {
    					IAtomContainer test = radii[r];

    						
    					/*for( IBond b : test.bonds() ) {
    						List<Integer> bondOrigins = (List<Integer>) b.getProperty( CDKSMARTSHyperstructureFitness.bondMolOriginType );
    					}*/
    					
    					boolean isGhost = ConvenienceTools.isGhostSubstructure(test);

    					
    					
    					int identifier = (Integer) test.getProperty( ConvenienceTools.origBondIndicesProperty );
    					
    					if( this.detailedGhostInfo && isGhost && ! uniqueSubstructures.containsKey(identifier) )
    						if( test instanceof IQueryAtomContainer ) {
    							QueryAtomContainer testQ = (QueryAtomContainer) test;
    							mcsGhostSMARTSSet.add( this.smaH.toSmarts( testQ ) );
    						} else {
    							try {
									mcsGhostSMARTSSet.add( this.sGenerator.create( test ) );
								} catch (CDKException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
    						}
    					
    					uniqueSubstructures.put( identifier, isGhost );
    					
    					
    				}
    			}
    			
    			
    			for( boolean b : uniqueSubstructures.values() ) {
    				if(b)
    					mcsGhostCount++;
    			}
    			
    			uniqueSubstructureCount = uniqueSubstructures.size();
    		}
		
		} else {
			bondMaps = new ArrayList<Map<IBond, IBond>>();
			atomMaps = new ArrayList<AtomAtomMapping>();
			atomIndexMaps = new ArrayList<Map<Integer, Integer>>();
			bondIndexMaps =  new ArrayList<List<int[]>>();
			
			fragmentSizes = new int[]{0};
		}
	}
	
	
	
	public double getSimilarityValue( SimilarityType type ) {
		
		if( type == SimilarityType.Tanimoto )
			return tanimoto;
		
		if( type == SimilarityType.Tversky )
			return tversky;
		
		
		return bondMaps.get(0).size();
	}
	
	
	
	//ExtendedIsomorphism comparator;
	private MCSMethods mcsMapper;
	private boolean useWeights;
	private boolean ghostSubstructures;
	//private boolean detailedGhostInfo;
	private int ghostRadius = 4;
	
	private ExtendedAlgorithm mappingAlgorithmName;
	private boolean useRaymondHeuristics;
	private boolean useRingHeuristics;
	private int topologicalDistanceLimit;
	private int expansionTimeOut;

	private SmartsHelper smaH;
	private SmilesGenerator sGenerator;
	private boolean detailedGhostInfo;
	
	public double tversky = -1.0;
	public double tanimoto = -1.0;
	
	
	public List<Map<IBond, IBond>> bondMaps;
	public List<List<int[]>> bondIndexMaps;
	//public List<Map<IAtom, IAtom>> atomMaps;
	public List<AtomAtomMapping> atomMaps;
	public List<Map<Integer, Integer>> atomIndexMaps;
	//public int refBondWeightSum;
	public int[] commonBondWeights = null;
	public String[] commonBondTopologies = null;
	//public int[] uniqueHSWeights = null;
	//public int[] uniqueDbWeights = null;
	public int[] fragmentSizes = null;
	public String mcsSMARTS = null;
	int uniqueSubstructureCount;
	int mcsGhostCount;
	public List<String> mcsGhostSMARTSSet; 
	public int modProdNodeCount;
	public long modProdConstructionTime;
	public long mcsExecTime;
	public double modProdEdgeDensity;
	
	
}