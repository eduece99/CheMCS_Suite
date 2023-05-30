package org.cisrg.mapping;


import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
//import org.knime.cisrg.hyperstructures.GAPlugins;
import org.cisrg.mapping.GenerateCompatibilityGraphEdges;
//import org.knime.cisrg.hyperstructures.GenerateCompatibilityGraphEdgesEdgeTypes;
//import org.knime.cisrg.hyperstructures.GenerateCompatibilityGraphMatrix;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
 
 



/**
 * 
 * @author edmund duesbury
 * @date 15th January 2015
 */
public abstract class CliqueDetection extends MCSMethods {

	/*
	 * Differences between this implementation and John Raymond's paper:
	 * 
	 * - MSI isn't used - we don't care if two molecules are very different, we still want the MCES.  Lower bound is therefore maximum clique size
	 * - We don't use the colouring algorithm, don't see the need if the partitioning is better
	 * 
	 * 
	 */
	
	 
	
	protected CliqueDetection( ModularProductOptions opts ) {
		super();
		
		this.mpOpts = opts;
	}
	
	
	protected CliqueDetection( ModularProductOptions opts, boolean swapTest ) {
		super();
		
		this.mpOpts = opts;
		this.swapTest = swapTest;
	}
	
	
	
	protected void saveClique( Collection<Integer> currentClique ) {
		
		//boolean deltaY = deltaYExchangeOccured( currentClique );
		
		// size check
		if( bestCliques.size() > 10000 )
			return;
		
		if( deltaYPossible ) {

			if( bestCliques.isEmpty() ) {
				if( ! ConvenienceTools.deltaYExchangeOccured( modProd, hsMol, qMol, currentClique ) ) {
					bestCliques.add( new ArrayList<Integer>(currentClique) );
					lowerBound = currentClique.size();
				}
			} else if( currentClique.size() >= bestCliques.get(0).size() ) {
				if( ! ConvenienceTools.deltaYExchangeOccured( modProd, hsMol, qMol, currentClique ) ) {
					if( currentClique.size() > bestCliques.get(0).size() ) {  // delete smaller cliques than newly-identified best
						bestCliques.clear();
						lowerBound = currentClique.size();
					}

					bestCliques.add( new ArrayList<Integer>(currentClique) );
				}
			}

		} else {

			if( bestCliques.isEmpty() ) {
				bestCliques.add( new ArrayList<Integer>(currentClique) );
				lowerBound = currentClique.size();
			} else if( currentClique.size() >= bestCliques.get(0).size() ) {

				if( currentClique.size() > bestCliques.get(0).size() ) {  // delete smaller cliques than newly-identified best
					bestCliques.clear();
					lowerBound = currentClique.size();
				}

				bestCliques.add( new ArrayList<Integer>(currentClique) );

			}

		}
		

		if( verbose ) {
		
			// debug stuff
			Map<IBond, IBond> commonSubgraphMap = new HashMap<IBond, IBond>();
			
			for( Integer mpn : bestCliques.get( bestCliques.size() - 1 ) ) {
				int[] node = modProd.getNodes().get(mpn);
				
				IBond hsBond = hsMol.getBond(node[0]);
				IBond qBond = qMol.getBond(node[1]);
	
				commonSubgraphMap.put(hsBond, qBond);
			}
			
			
			IAtomContainer sg = ConvenienceTools.createCommonSubgraph(hsMol, qMol, commonSubgraphMap);
			try {
				System.out.println( "current clique - " + sGenerator.create(sg) );
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	

	protected abstract void findCliques();
	
	

	public void search( IAtomContainer graph1, IAtomContainer graph2 ) {
		
		hsMol = graph1;
		qMol = graph2;
		
		boolean swapMolecules = false;
		if( swapTest ) {
			if( hsMol.getBondCount() > qMol.getBondCount() ) {
				swapMolecules = true;
			}
		}
		
		if( swapMolecules ) {
 				IAtomContainer tempMol = qMol;
 				qMol = hsMol;
 				hsMol = tempMol;
		}
		
		// Here we test for the presence of "triangles" to account for delta-Y exchanges
		// that is, as they mess up the translation from max clique in modular product, to the MCES
		//long timeThing = System.currentTimeMillis();
		//ConvenienceTools.countRings(hsMol);
		try {
			
			IAtomContainer cyclopropane = sParser.parseSmiles("*1**1");
			

			boolean cp1 = (hsMol instanceof IQueryAtomContainer) ? true : ConvenienceTools.isSubgraph(cyclopropane, hsMol);
			boolean cp2 = ConvenienceTools.isSubgraph(cyclopropane, qMol);
			
			//if( uit.isSubgraph(hsMol, cyclopropane) || uit.isSubgraph(qMol, cyclopropane) ) {
			if( cp1 || cp2 ) {
				deltaYPossible = true;
			} else {
				deltaYPossible = false;
			}
			
			//timeThing = System.currentTimeMillis() - timeThing;
			//System.out.println("time taken test - " + timeThing + " " + deltaYPossible );
		} catch (InvalidSmilesException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} 
		//algTime = System.currentTimeMillis();
		modProdTime = System.currentTimeMillis();
		
		modProd = null;
		try {
			if( isConnected ) {
				modProd = new GenerateCompatibilityGraphEdgesEdgeTypes(
						hsMol, 
						qMol, 
						matchBonds, 
						mpOpts.useRaymondHeuristics, 
						mpOpts.useRingHeuristics, 
						mpOpts.topoDistanceLimit, 
						null, 
						true 
				);
			} else {
				modProd = new GenerateCompatibilityGraphEdges(
						hsMol, 
						qMol, 
						matchBonds, 
						mpOpts.useRaymondHeuristics, 
						mpOpts.useRingHeuristics, 
						mpOpts.topoDistanceLimit, 
						null, 
						true 
				);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		modProdTime = System.currentTimeMillis() - modProdTime;
		int modProdNodeCount = modProd.getNodes().size();
		modProdEdgeDensity = (modProd.numberOfEdges ) / (double)( modProdNodeCount * (modProdNodeCount - 1) );
		

		lowerBound = 0;  // initialise lower bound
		
		bestCliques = new LinkedList<List<Integer>>();
		
		
		
		
		// null solution place-holder
		if( modProd.getNodes().size() > 0 && modProd.numberOfEdges > 0 ) {
			findCliques( );
		} else {
			bestCliques.add( new ArrayList<Integer>() );
		}
		
		
		// re-swap here
		
		if( swapMolecules ) {
			IAtomContainer tempMol = qMol;
			qMol = hsMol;
			hsMol = tempMol;
		}
					
					
		// cliques to solutions
		for( List<Integer> clique : bestCliques ) {
			
			List<int[]> pairList = new ArrayList<int[]>( clique.size() );
			for( Integer index : clique ) {
				int[] pair = modProd.getNodes().get(index);
				
				int[] newPair = new int[2];
				
				if( swapMolecules ) {
					newPair[0] = pair[1];
					newPair[1] = pair[0];
				} else {
					newPair[0] = pair[0];
					newPair[1] = pair[1];
				}
				
				pairList.add(newPair);
			}
			mcsBondIndexIsomorphisms.add( pairList );
			
			Map<IBond, IBond> commonSubgraphMap = new HashMap<IBond, IBond>();
			
			
			for( Integer mpn : clique ) {
				int[] node = modProd.getNodes().get(mpn);
				
				IBond hsBond;
				IBond qBond;
				
				if( swapMolecules ) {
					 hsBond = hsMol.getBond(node[1]);
					 qBond = qMol.getBond(node[0]);
				} else {
					 hsBond = hsMol.getBond(node[0]);
					 qBond = qMol.getBond(node[1]);
				}
				
				
				/*System.out.println( 
						hsMol.getAtomNumber( hsBond.getAtom(0) ) + "," +
						hsMol.getAtomNumber( hsBond.getAtom(1) )  + "=" +
						qMol.getAtomNumber( qBond.getAtom(0) ) + "," +
						qMol.getAtomNumber( qBond.getAtom(1) )
				);*/
				
				commonSubgraphMap.put(hsBond, qBond);
			}
			
			
			
			
			
			Map<IAtom, IAtom> commonAtomMap = createAtomMapFromBondMap(hsMol, qMol, commonSubgraphMap);
			
			// query to hyperstructure mapping
			List<Integer> mapping = atomMapToChromosome(hsMol, qMol, commonAtomMap);
			//System.out.println( "mapping - " + commonAtomMap  );
			//System.out.println( "chromosome - " + mapping  );
			
			mcsAtomIsomorphisms.add( mapping );
			mcsBondIsomorphisms.add( commonSubgraphMap );
			
			
			//algTime = System.currentTimeMillis() - algTime;
		}
		
		
		
		
		
		if( verbose ) {
			System.out.println( "\n max clique size = " + lowerBound + " " + 1L );
			//System.out.println( "upperbound = " + upperBound );
			
			System.out.println( "Modular Product nodes = " + modProdNodeCount );
			System.out.println( "Modular Product edges = " + modProd.numberOfEdges );
			System.out.println( "Modular Product edge density = " + modProdEdgeDensity );
			System.out.println( "Modular Product creation time = " + modProdTime );
		}
		//System.out.println( "buildup time - " + buildupTime );
	}
	
	
	
	protected int getFirstCliqueSize() {
		return bestCliques.get(0).size();
	}
	
	
	
	
	
	// descending order of cardinality
	protected Comparator<List> sizeDescendingComparator = new Comparator<List>() {
			 public int compare(List o1, List o2) {

		            int dComp;
		            // sort in descending order of degree first
		            if( o1.size() == o2.size() ) {
		            	dComp = 0;
		            } else if( o1.size() > o2.size() ) {
		            	dComp = -1;
		            	//Collections.swap( sortedVertices2, o1, o2 );
		            } else {
		            	dComp = 1;
		            	//Collections.swap( sortedVertices2, o2, o1 );
		            }
		            
		            return dComp;
			 } 
		};
		
			

		
	protected int numberOfSteps = 0;
	//protected long startTime;
	protected boolean timeOut = false;	
	protected boolean isConnected = false;
	protected boolean swapTest = false;
	
	protected List< List<Integer> > bestCliques;
	protected int lowerBound;
	
	
	
	protected double MinSimilarityIndex = 0.01;  // currently unused
	protected boolean deltaYPossible;
	protected ModularProductOptions mpOpts;



	public static class ModularProductOptions {
		
		public ModularProductOptions( boolean raymondHeuristics, boolean ringHeuristics, int topoDistLimit ) {
			this.useRaymondHeuristics = raymondHeuristics;
			this.useRingHeuristics = ringHeuristics;
			this.topoDistanceLimit = topoDistLimit;
		}
		
		
		protected boolean useRaymondHeuristics = false;
		protected boolean useRingHeuristics = false;
		protected int topoDistanceLimit = -1;
		
	}

	
	
	
}
