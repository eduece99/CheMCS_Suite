package org.cisrg.mapping;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
//import org.knime.cisrg.hyperstructures.GAPlugins;
import org.openscience.cdk.exception.CDKException;
 



public class CarraghanPardalosMCS extends CliqueDetection {

	
	public CarraghanPardalosMCS(ModularProductOptions opts) {
		super( opts );
	}
	
	protected void orderNodes() {
		
		// variable initialisation
		nNodes = modProd.getNodes().size();
		mat = new ArrayList<List<Boolean>>( nNodes );
		int[] edge = new int[ nNodes ];  // degree values for nodes
		actNode = new ArrayList<Integer>( nNodes );
		
		// Convert Modular Product into adjacency matrix
		List<Collection<Integer>> adjList = modProd.getAdjacencyList();
		
		
		for( int r = 0; r < nNodes; r++ ) {
			List<Boolean> amRow = new ArrayList<Boolean>( nNodes );
			
			// initial assignments are 'off'
			for( int c = 0; c < nNodes; c++ ) {
				amRow.add( false );
			}
			
			// now to populate with adjacency information
			for( Integer index : adjList.get(r) ) {
				amRow.set( index, true );
			}
			
			mat.add( amRow );
			
			// maintain pointers to original matrix
			actNode.add( r );
		}
		//actNode.add( nNodes );
		
		double density = (modProd.numberOfEdges ) / (double)( nNodes * (nNodes - 1) );
		
		// perform ordering (if needed)
		if( density >= 0.4 ) {
			
			// count degrees of nodes
			for( int row = 0; row < nNodes; row++ ) {
				edge[row] = 0;
				
				for( int col = 0; col < nNodes; col++ ) {
					if( mat.get(row).get(col) ) {
						edge[row]++;  
					}
				}
			}
			
			
			int minNode = 0;
			
			for( int node = 0; node < nNodes - 2; node++ ) {
				int min = nNodes;
				
				// find the node of smallest degree in the current adjacency matrix
				for( int row = node; row < nNodes; row++ ) {
					if( edge[row] < min ) {
						min = edge[row];
						minNode = row;
					}
				}
				
				// re-assign degree
				edge[minNode] = edge[node];
				
				// swap node pointers
				if( node != minNode ) {
					int temp = actNode.get(node);
					actNode.set( node, actNode.get(minNode) );
					actNode.set( minNode, temp );
				
					// adjust adjacency matrix adjacency values
					for( int row = 0; row < nNodes; row++ ) {
						boolean tempLog = mat.get(row).get(node);
						mat.get(row).set(node, mat.get(row).get(minNode));
						mat.get(row).set(minNode, tempLog);
					}
					
					for( int col = 0; col < nNodes; col++ ) {
						boolean tempLog = mat.get(node).get(col);
						mat.get(node).set(col, mat.get(minNode).get(col) );
						mat.get(minNode).set(col, tempLog);
					}
				}
				
				// update degree values (removal of said node basically)
				for( int col = node; col < nNodes; col++ ) {
					if( mat.get(node).get(col) ) {
						edge[col]--;
					}
				}
			}
			
		}
		
	}
	
	
	protected void findCliquesCP() throws CDKException {
		
		int[] start = new int[ nNodes ];  // node currently being expanded, at depth d
		int[] last = new int[ nNodes ];		// last node to (possibly) be expanded at depth d
		List<Integer> mClique = new ArrayList<Integer>( nNodes );
		int d = 1;  // depth
		int dTemp;
		
		// initialisation of variables
		start[1] = 0;
		last[1] = nNodes;
		int[][] adj = new int[ nNodes + 1 ][ nNodes + 1 ];
		for( int col = 0; col < nNodes; col++ ) {
			adj[1][col] = col;
		}
		
		while( d > 0 ) {
			
			// XXX  time-out condition
			++numberOfSteps;
			
			if( numberOfSteps % 10000 == 0 ) {
				long currentTime = System.currentTimeMillis();
				if( currentTime - mcsStartTime > expansionTimeLimitMs ) {
					throw new CDKException("Error - time limit of " + expansionTimeLimitMs + " ms reached!");
				}
			}
			
			start[d]++;
			if( (d + last[d] - start[d]) > lowerBound ) {
				dTemp = d;
				d++;
				start[d] = 0;
				last[d] = 0;
				
				// determine the node to use at the next depth
				for( int col = start[dTemp]; col <= last[dTemp]; col++ ) {
					if( mat.get( adj[dTemp][start[dTemp]] ).get( adj[dTemp][col]) ) {
						last[d]++;
						adj[d][ last[d] ] = adj[dTemp][col];
					}
				}
				
					/*
					 * if the next depth doesn't contain any nodes, see if a new
					 * maximum clique has been found, and return to the previous depth
					 */
				if( last[d] == 0 ) {
					d--;
						
					if( d > lowerBound ) {
						//maxCliqueSize = d;
						mClique.clear();
						
						for( int col = 1; col <= d; col++ ) {
							mClique.add( actNode.get( adj[col][ start[col] ] ) );
						}
						
						saveClique(mClique);
					}
				}
			} else { // prune, further expansion would not find a better incumbent
				d--;
			}
		}
		
		System.out.println( "max clique size " + lowerBound );
		//bestCliques.add( mClique );
	}

	@Override
	protected void findCliques() {
		orderNodes();
		
		try {
			findCliquesCP();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	
	int nNodes;
	List<Integer> actNode;		// "actual" nodes - translation from modded adj matrix to original graph node indices
	List<List<Boolean>> mat;	// adjacency matrix
	//int maxCliqueSize = 0;
	
}
