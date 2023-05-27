package org.cisrg.mapping;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
//import org.knime.cisrg.hyperstructures.GAPlugins;
import org.openscience.cdk.exception.CDKException;






public class BronKerboschCazals extends CliqueDetection {

	
	
	public BronKerboschCazals(ModularProductOptions opts) {
		super( opts );
	}
	
	
	/**
	 * 
	 * 
	 * @param r Current Clique set
	 * @param p Pivot vertices
	 * @param x Exclusion set
	 * @throws CDKException
	 */
	private void findCliquesBK( Collection<Integer> r, Collection<Integer> p, Collection<Integer> x ) throws CDKException {
		
		// time-out condition
		++numberOfSteps;
		
		if( numberOfSteps % 1000 == 0 ) {
			long currentTime = System.currentTimeMillis();
			if( currentTime - mcsStartTime > expansionTimeLimitMs ) {
				throw new CDKException("Error - time limit of " + expansionTimeLimitMs + " ms reached!");
			}
		}
		
		if( p.isEmpty() && x.isEmpty() ) {
			
			saveClique(r);

			//System.out.println( r );
		} else {
			
			Collection<Integer> pUx = ConvenienceTools.createCollection( p );
			pUx.addAll( x );
			Integer pivot = pUx.iterator().next();
			for( Integer pV : pUx ) {
				// choose the pivot vertex - maximum degree of potential pivots
				if( modProd.getAdjacencyList().get(pV).size() > modProd.getAdjacencyList().get(pivot).size() ) {
					pivot = pV;
				}
			}
			
			
			//Collection<Integer> pivotNonNeighbours = new HashSet<Integer>( pUx );
			Collection<Integer> pivotNonNeighbours = pUx;
			
			if( pivot != null ) {
				Collection<Integer> pivotNeighbours = modProd.intersectNodeNeighbours( pivot, p );
				//for( Integer i : pivotNeighbours ) { pivotNonNeighbours.remove( i ); } 
				pivotNonNeighbours.removeAll(pivotNeighbours);
				//pivotNonNeighbours = new HashSet<Integer>( pivotNonNeighbours );
			}
			
			// look at vertices not adjacent to the pivot
			for( Integer v : pivotNonNeighbours ) {
				
				//p.remove(v);
				
				Collection<Integer> rNew = new HashSet<Integer>( r );
				rNew.add( v );
				
				// if delta-Y exchange occurs, save the non-exchanged clique
				if( deltaYPossible && ConvenienceTools.deltaYExchangeOccured( modProd, hsMol, qMol, new ArrayList<Integer>(rNew) ) ) {
					//bestCliques.add( new ArrayList<Integer>( r ) );
					saveClique( r );
					continue;
				}
				
				/*
				System.out.println( "clique size - " + rNew.size() );
				System.out.println( "p size - " + p.size() );*/
				
				Collection<Integer> pNew = ConvenienceTools.createCollection( p.size() + 1 );
				Collection<Integer> pIntersected = modProd.intersectNodeNeighbours(v, p );
				//pNew.clear();
				pNew.addAll( pIntersected );
				//for( Integer vP : pIntersected ) { pNew.add( vP ); }
				
				Collection<Integer> xNew = ConvenienceTools.createCollection( modProd.getNodes().size() );
				Collection<Integer> xIntersected = modProd.intersectNodeNeighbours(v, x );
				//xNew.clear();
				xNew.addAll( xIntersected );
				//for( Integer vX : xIntersected ) { xNew.add( vX ); }
				
				findCliquesBK( rNew, pNew, xNew );
				
				x.add( v );
			}
			
		}
		
	}
	
	
	@Override
	protected void findCliques() {
		
		Collection<Integer> initialNodes = ConvenienceTools.createCollection( modProd.getNodes().size() );
		for( int n = 0; n < modProd.getNodes().size(); n++ ) { initialNodes.add(n); }
		
		try {
			findCliquesBK( 
					ConvenienceTools.createCollection( modProd.getNodes().size() ), 
					initialNodes, 
					ConvenienceTools.createCollection( modProd.getNodes().size() ) 
			);
		} catch (CDKException e) {
			e.printStackTrace();
		}
	}
	




	
}
