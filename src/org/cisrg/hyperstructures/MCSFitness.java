package org.cisrg.hyperstructures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.cisrg.knime.ExtendedIsomorphism;
import org.cisrg.mapping.ConvenienceTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.smsd.Isomorphism;

/**
 * Uses the CDK API to convert a given map to a bond map, and uses the number of correctly-assigned bonds as the fitness measure
 * 
 * @author edmund
 *
 */
public class MCSFitness extends GAPlugins {

	public MCSFitness( IAtomContainer hsMol, IAtomContainer qMol ) {
		this.qMol = hsMol;
		this.hsMol = qMol;
	}
	
	@Override
	public double fitness(ArrayList<Integer> mapChr) {
		
	HashMap<IAtom, IAtom> atomMap = new HashMap<IAtom, IAtom>();
 	   
 	   //System.out.println("stuff: " + mapChr.size() + " " + mapper.getQueryMol().getAtomCount() + " " + mapper.getMainMol().getAtomCount() + " | " + mapChr );
 	   
 	   // translate mapping from array to hash
 	   for( int m = 0; m < mapChr.size(); m++ ) {
 		   if( mapChr.get(m) >= 0 ) {
 			   if( hsMol.getAtom( m ) != null && qMol.getAtom( mapChr.get(m) ) != null )
 				   atomMap.put( hsMol.getAtom( m ), qMol.getAtom(  mapChr.get(m) ) );
 		   }
 	   }
		 
 	  Map<IBond, IBond> bondMap = ConvenienceTools.makeBondMapOfAtomMap( hsMol, qMol, atomMap );
 	 //System.out.println( bondMap.size() + " " + bondMap );
 	  return bondMap.size() / (double) (qMol.getBondCount() + 0.00001) ;
	}

	private IAtomContainer hsMol, qMol;
	
}
