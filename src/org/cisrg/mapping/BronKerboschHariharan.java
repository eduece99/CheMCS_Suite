package org.cisrg.mapping;

import java.util.ArrayList;
import java.util.Collection;
//import org.knime.cisrg.hyperstructures.BitSetExtended;
//import org.knime.cisrg.hyperstructures.GAPlugins;
import org.cisrg.mapping.GenerateCompatibilityGraphEdges;
import org.cisrg.mapping.GenerateCompatibilityGraphEdgesEdgeTypes;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
 




/**
 * Modification of the Bron-Kerbosch algorithm to find maximum c-cliques.  Original idea by Hariharan et al (2011) - MultiMCS
 * 
 * Currently finds maximum c-cliques.  Can be easily modified to find maximal ones instead (saveClique function)
 * 
 * @author edmund duesbury
 *
 */
public class BronKerboschHariharan extends CliqueDetection {

	
	GenerateCompatibilityGraphEdgesEdgeTypes modProdC = null;
	
	
	public BronKerboschHariharan(ModularProductOptions opts) {
		
		super( opts );
		
		isConnected = true;  // needed for modular product set-up (additional edge info)
	}
	
	
	/**
	 * 
	 * @param r
	 * @param n - set of vertices
	 * @param x - exclusion set
	 * @param ng - green neighbours of clique
	 * @param nr - red neighbours of clique
	 * 
	 * @throws CDKException
	 */
	private void findCliquesBK( Collection<Integer> r, Collection<Integer> n, Collection<Integer> ng, Collection<Integer> nr, Collection<Integer> x ) throws CDKException {
		
		// time-out condition
		++numberOfSteps;
		
		if( numberOfSteps % 1000 == 0 ) {
			long currentTime = System.currentTimeMillis();
			if( currentTime - mcsStartTime > expansionTimeLimitMs ) {
				throw new CDKException("Error - time limit of " + expansionTimeLimitMs + " ms reached!");
			}
		}
		
		if( n.isEmpty() ) {
			
			saveClique(r);

			//System.out.println( r );
		} else {
			
			
			
			 
			for( Integer v : n ) {
				
				// avoid exclusion set
				if( x.contains(v) ) {
					continue;
				}
				
				
				Collection<Integer> rNew = ConvenienceTools.createCollection( r );
				rNew.add( v );
				
				
				
				// if delta-Y exchange occurs, save the non-exchanged clique
				if( deltaYPossible && ConvenienceTools.deltaYExchangeOccured( modProdC, hsMol, qMol, new ArrayList<Integer>(rNew) ) ) {
					//bestCliques.add( new ArrayList<Integer>( r ) );
					saveClique( r );
					continue;
				}
				
				// green edges (from vertex)
				Collection<Integer> greenNeighbours = modProdC.getCEdges().get(v);
				Collection<Integer> ngNew = ConvenienceTools.createCollection(ng);
				ngNew.addAll( greenNeighbours );
				ngNew.removeAll( rNew );  // remove current clique nodes from expansion candidates
				
				Collection<Integer> redNeighbours = modProdC.getREdges().get(v);
				Collection<Integer> nrNew = ConvenienceTools.createCollection(nr);
				nrNew.addAll( redNeighbours );
				
				
				Collection<Integer> ngMinusNr = ConvenienceTools.createCollection(ngNew);
				ngMinusNr.removeAll(nrNew);
				
				
				Collection<Integer> xNew = ConvenienceTools.createCollection( x );
				
				
				findCliquesBK( rNew, ngMinusNr, ngNew, nrNew, xNew );
				
				x.add( v );
			}
			
		}
		
	}
	
	
	@Override
	protected void findCliques() {
		
		modProdC = ( GenerateCompatibilityGraphEdgesEdgeTypes ) modProd;
		
		Collection<Integer> initialNodes = ConvenienceTools.createCollection( modProdC.getNodes().size() );
		for( int n = 0; n < modProdC.getNodes().size(); n++ ) { initialNodes.add(n); }
		
		
		try {
			findCliquesBK( 
					ConvenienceTools.createCollection( modProdC.getNodes().size() ), 
					initialNodes, 
					ConvenienceTools.createCollection( modProdC.getNodes().size() ), 
					ConvenienceTools.createCollection( modProdC.getNodes().size() ), 
					ConvenienceTools.createCollection( modProdC.getNodes().size() ) 
			);
		} catch (CDKException e) {
			e.printStackTrace();
		}
	}
	




	
}
