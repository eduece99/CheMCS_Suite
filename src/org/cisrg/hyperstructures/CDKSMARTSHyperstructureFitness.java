package org.cisrg.hyperstructures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.cisrg.ambit.AliphaticSymbolQueryAtom;
import org.cisrg.ambit.AromaticSymbolQueryAtom;
import org.cisrg.ambit.RingClosure;
import org.cisrg.ambit.SingleOrAromaticBond;
import org.cisrg.ambit.SmartsAtomExpression;
import org.cisrg.ambit.SmartsBondExpression;
import org.cisrg.knime.ExtendedIsomorphism;
import org.cisrg.mapping.ConvenienceTools;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.isomorphism.matchers.OrderQueryBond;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.smarts.AliphaticSymbolAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AnyAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AnyOrderQueryBond;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticQueryBond;
import org.openscience.cdk.isomorphism.matchers.smarts.AromaticSymbolAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.AtomicNumberAtom;
import org.openscience.cdk.isomorphism.matchers.smarts.SMARTSBond;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

//import ambit2.smarts.AromaticSymbolQueryAtom;
import ambit2.smarts.SingleBondAromaticityNotSpecified;
import ambit2.smarts.SmartsConst;
import ambit2.smarts.SmartsExpressionToken;

public class CDKSMARTSHyperstructureFitness extends CDKHyperstructureFitness {

	public CDKSMARTSHyperstructureFitness( IQueryAtomContainer hs, IAtomContainer qmol ) {
		super();
		
		this.hs = hs;
		this.queryMol = qmol;
		
		/*if( bondFrequencies ) {
			for( IBond b : hs.bonds() ) {
				if( b.getProperty(bondFrequencyType) == null )
					b.setProperty(bondFrequencyType, 1);
				
				if( ConvenienceTools.isRingBond( b ) ) {
					b.setProperty(topologyType, "r");
				} else {
					b.setProperty(topologyType, "c");
				}
				
				if( b.getProperty(bondMolOriginType) == null ) {
					List<Integer> initial = new ArrayList<Integer>( 10 );
					initial.add(0);
					
					b.setProperty(bondMolOriginType, initial );  	
				}
			}
		}*/
		
		//configureBondOrders(hs);
	}
	
	public CDKSMARTSHyperstructureFitness( IQueryAtomContainer hs, IAtomContainer qmol, boolean ringEnforcement, boolean bondFrequencies, int qMolId ) {
		super();
		
		this.hs = hs;
		this.queryMol = qmol;
		
		this.ringEnforcement = ringEnforcement;
		this.bondFrequencies = bondFrequencies;
		this.queryMolId = qMolId;
		
		ConvenienceTools.countRings(hs);
		
		if( bondFrequencies ) {
			for( IBond b : hs.bonds() ) {
				if( b.getProperty(bondFrequencyType) == null )
					b.setProperty(bondFrequencyType, 1);
				
					if( ConvenienceTools.isRingBond( b ) ) {
						b.setProperty(topologyType, "r");
					} else {
						b.setProperty(topologyType, "c");
					}
					
				if( b.getProperty(bondMolOriginType) == null ) {
					List<Integer> initial = new ArrayList<Integer>( 10 );
					initial.add( 0 );
					
					b.setProperty(bondMolOriginType, initial );  	
				}
			}
		}
		
		//configureBondOrders(hs);
	}
	
	/**
	 *  This mainly concerns "hacking" the hyperstructure to be exportable to the MDLV2000 format
	 *  
	 *  Thus, giving atoms dummy symbols and bonds dummy orders, is how this works.
	 *  
	 * @param qac
	 */
	public static void configureBondOrders( IQueryAtomContainer qac ) {
		for( IBond b : qac.bonds() ) {
			if( b.getOrder() == null || b.getOrder() == Order.UNSET ) {
				b.setOrder( Order.SINGLE );
			}
		}
		
		for( IAtom a : qac.atoms() ) {
			if( a.getSymbol() == null ) {
				a.setSymbol( "Si" );
			}
		}
		/*
		// remove hydrogens as they cause errors in CDK 1.5.2 upon subgraph matching
		 * No longer needed, using CDK 1.5.4
		for( IAtom at : qac.atoms() ) {
			if( at.getSymbol() != null && at.getSymbol().equals("H") ) {
				qac = new QueryAtomContainer( AtomContainerManipulator.removeHydrogens( qac ), qac.getBuilder() );
				break;
			}
		}
		*/
	}
	
	private int bondOrderToAMBITBond( Order order ) {
		switch( order ) {
		case SINGLE:
			return SmartsConst.BT_SINGLE;
		case DOUBLE:
			return SmartsConst.BT_DOUBLE;
		case TRIPLE:
			return SmartsConst.BT_TRIPLE;
		default:
			return SmartsConst.BT_AROMATIC;
		}
	}
	
	public void setRingBonds( ) {
		
		
		IQueryAtomContainer nHs = (IQueryAtomContainer) hs;
		

		IBond[] newBonds = new IBond[ nHs.getBondCount() ] ;
		
		for( int b = 0; b < nHs.getBondCount(); b++ ) {
			IBond oldBond = nHs.getBond(b);
			
			IBond newBond = oldBond;
			
			IAtom atom1 = oldBond.getAtom(0);
			IAtom atom2 = oldBond.getAtom(1);
	
			
			/*if( oldBond instanceof RingBond ) {
				newBond = new RingBond();
			} else*/ 
			
			//boolean isAromatic = ConvenienceTools.isAromatic( oldBond );
			
			//newBond = createSMARTSRingChainBond(atom1, atom2, isRingBond( oldBond ), isAromatic, oldBond.getOrder() );
			newBond = createSMARTSRingChainBond(atom1, atom2, oldBond, oldBond );
			newBond.setProperties( oldBond.getProperties() );
			//System.out.println( newBond.getProperties() );
		
			newBond.setAtom(atom1, 0);
			newBond.setAtom(atom2, 1);
			
			newBond.setProperty(origBondType, oldBond);
			
			newBonds[b] = newBond;
		}

		
		nHs.setBonds( newBonds );
		
		
	}
	
	public IQueryAtomContainer getAllBondMatchHyperstructure( boolean ignoreLabels ) {
		
		IQueryAtomContainer nHs = (IQueryAtomContainer) hs;
	

		IBond[] newBonds = new IBond[ nHs.getBondCount() ] ;
		IAtom[] newAtoms = new IAtom[ nHs.getAtomCount() ] ;
		
		for( int b = 0; b < nHs.getBondCount(); b++ ) {
			IBond oldBond = nHs.getBond(b);
			
			IBond newBond = oldBond;
			
			IAtom atom1 = oldBond.getAtom(0);
			IAtom atom2 = oldBond.getAtom(1);
			
			int index1 = nHs.getAtomNumber(oldBond.getAtom(0));
			int index2 = nHs.getAtomNumber(oldBond.getAtom(1));
			
			if( ignoreLabels ) {
				
				
				if( newAtoms[index1] == null ) {
					atom1 = new AnyAtom( hs.getBuilder() );
				} else {
					atom1 = newAtoms[index1];
				}
				
				if( newAtoms[index2] == null ) {
					atom2 = new AnyAtom( hs.getBuilder() );
				} else {
					atom2 = newAtoms[index2];
				}
				
				atom1.setSymbol("C");
				atom2.setSymbol("C");
				
				newAtoms[ nHs.getAtomNumber(oldBond.getAtom(0)) ] = atom1;
				newAtoms[ nHs.getAtomNumber(oldBond.getAtom(1)) ] = atom2;
			} else {
				if( newAtoms[index1] == null && ( atom1 instanceof AromaticSymbolQueryAtom || atom1 instanceof AliphaticSymbolQueryAtom ) ) {
					String symbol1 = atom1.getSymbol();
					atom1 = new AtomicNumberAtom( ConvenienceTools.atomSymbolToNumber( atom1 ), DefaultChemObjectBuilder.getInstance() );
					atom1.setSymbol( symbol1 );
				} else if( newAtoms[index1] != null ) {
					atom1 = newAtoms[index1];
				} else {
					atom1 = nHs.getAtom(index1);
				}
					
				if( newAtoms[index2] == null && ( atom2 instanceof AromaticSymbolQueryAtom || atom2 instanceof AliphaticSymbolQueryAtom ) ) {
					String symbol2 = atom2.getSymbol();
					atom2 = new AtomicNumberAtom( ConvenienceTools.atomSymbolToNumber( atom2 ), DefaultChemObjectBuilder.getInstance() );
					atom2.setSymbol( symbol2 );
				} else if( newAtoms[index2] != null ) {
					atom2 = newAtoms[index2];
				} else {
					atom2 = nHs.getAtom(index2);
				}
				
				newAtoms[ index1 ] = atom1;
				newAtoms[ index2 ] = atom2;
			}
				//nHs.setAtom( nHs.getAtomNumber(oldBond.getAtom(0)), atom1 );
				//nHs.setAtom( nHs.getAtomNumber(oldBond.getAtom(1)), atom2 );
			
			
			
			/*if( oldBond instanceof RingBond ) {
				newBond = new RingBond();
			} else*/ 
			
			/*
			if( oldBond instanceof SingleOrAromaticBond ) {
				if( ConvenienceTools.isAromatic( oldBond.getAtom(0) )  && ConvenienceTools.isAromatic( oldBond.getAtom(1) ) ) {
					oldBond = new AromaticQueryBond( (IQueryAtom) oldBond.getAtom(0), (IQueryAtom) oldBond.getAtom(1), Order.SINGLE, DefaultChemObjectBuilder.getInstance() );
				}
			}
			*/
			/*
			if( newBond instanceof AromaticQueryBond || newBond instanceof AromaticOrSingleQueryBond || newBond instanceof SingleOrAromaticBond ) {
				//newBond = new AromaticQueryBond( (IQueryAtom) atom1, (IQueryAtom) atom2, Order.TRIPLE, hs.getBuilder() );
				//newBond.setOrder( Order.TRIPLE );
				atom1.setFlag(CDKConstants.ISAROMATIC, true);
				atom2.setFlag(CDKConstants.ISAROMATIC, true);
				//newBond.setProperty(manyOrderProp, 5);
			} 
			*/
			{
				newBond = new AnyOrderQueryBondNotNull( (IQueryAtom) atom1, (IQueryAtom) atom2, oldBond.getOrder(), hs.getBuilder() );
			}
			newBond.setAtom(atom1, 0);
			newBond.setAtom(atom2, 1);
			
			newBond.setProperty(origBondType, oldBond);
			newBond.setProperty(bondFrequencyType, oldBond.getProperty(bondFrequencyType) );
			newBond.setProperty(topologyType, oldBond.getProperty(topologyType) );
			newBond.setProperty(bondMolOriginType, oldBond.getProperty( bondMolOriginType ));
			//System.err.println( "oldbond = " + oldBond );
			newBonds[b] = newBond;
		}
		
		//if( ignoreLabels )
		nHs.setAtoms( newAtoms );
		
		nHs.setBonds( newBonds );
		
		
		
		return nHs;
	}
	
	/*
	public IQueryAtomContainer getAllAtomMatchHyperstructure() {
		
		IQueryAtomContainer nHs = (IQueryAtomContainer) hs;
		

		IAtom[] newBonds = new IAtom[ nHs.getAtomCount() ] ;
		
		for( int b = 0; b < nHs.getAtomCount(); b++ ) {
			IAtom oldBond = nHs.getAtom(b);
			IAtom newBond = new AnyAtom( hs.getBuilder() );
			
			newBonds[b] = newBond;
		}
		
		nHs.setAtoms( newBonds );
		
		return nHs;
	}
	*/
	
public IQueryAtomContainer restoreHyperstructureBonds( boolean ignoreLabels ) {
		
		IQueryAtomContainer nHs = (IQueryAtomContainer) hs;
		/*
		try {
			nHs = (IQueryAtomContainer) hs.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		*/
		IBond[] newBonds = new IBond[ nHs.getBondCount() ] ;
		IAtom[] newAtoms = new IAtom[ nHs.getAtomCount() ] ;
		HashMap<IAtom, IAtom> oldToNewAtoms = new HashMap<IAtom, IAtom>();
		
		for( int a = 0; a < nHs.getAtomCount(); a++ ) {
			IAtom oldAtom = nHs.getAtom(a);
			IAtom newAtom = null;
			
			if( ignoreLabels ) {
				newAtom = new AnyAtom( hs.getBuilder()) ;
				oldToNewAtoms.put( oldAtom, newAtom );
			} else  {
				if( oldAtom instanceof SmartsAtomExpression ) {
					newAtom = oldAtom;
				} else if( oldAtom.getProperty(atomicNumberProp) != null ) {
					//System.out.println( "new atom number symbol - " + oldAtom.getSymbol() );
					newAtom = new AtomicNumberAtom( ConvenienceTools.atomSymbolToNumber(oldAtom), hs.getBuilder() );
					newAtom.setSymbol( PeriodicTable.getSymbol( newAtom.getAtomicNumber() ) );
					oldToNewAtoms.put( oldAtom, newAtom );
					//newBond.setSymbol("C");
					//System.err.println( "bhah " + nHs.getAtomNumber(newBond) );
				} else if ( oldAtom.getProperty(aromaticAtomProp) != null ) {
					newAtom = new AromaticSymbolQueryAtom( hs.getBuilder() );
					//newAtom = new AromaticSymbolQueryAtom();
					newAtom.setSymbol( oldAtom.getSymbol() );
					oldToNewAtoms.put( oldAtom, newAtom );
				} else {
					newAtom = oldAtom;
				}
			}
			newAtoms[a] = newAtom;
		}
		
		nHs.setAtoms( newAtoms );
		
		
		for( int b = 0; b < nHs.getBondCount(); b++ ) {
			IBond oldBond = nHs.getBond(b);
			IBond newBond = oldBond;
			
			IAtom oldAtom1 = oldBond.getAtom(0);
			IAtom oldAtom2 = oldBond.getAtom(1);
			
			if( oldToNewAtoms.containsKey(oldAtom1) )
				oldAtom1 = oldToNewAtoms.get(oldAtom1);
			
			if( oldToNewAtoms.containsKey(oldAtom2) )
				oldAtom2 = oldToNewAtoms.get(oldAtom2);
			
			
			
			if( oldBond.getProperty(manyOrderProp) != null ) {
				System.err.println("GOATS");
				newBond = new AnyOrderQueryBondNotNull( (IQueryAtom) oldAtom1, (IQueryAtom) oldAtom2, oldBond.getOrder(), hs.getBuilder() );
				
				
				
				//if( oldBond.getProperty( topologyType ) != null )
					//newBond.setProperty(topologyType, oldBond.getProperty( topologyType ));
				
			} else if( ! ringEnforcement && oldBond.getProperty(origBondType) != null ) {
				//newBond = new OrderQueryBond( (IQueryAtom) oldAtom1, (IQueryAtom) oldAtom2, oldBond.getOrder(), hs.getBuilder() );
				
				try {
					IBond tempBond = oldBond.getProperty(origBondType);
					newBond = (IBond) tempBond.clone();
					newBond.setAtom( oldAtom1, 0);
					newBond.setAtom( oldAtom2, 1);
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				//newBond = new OrderQueryBond( (IQueryAtom) oldAtom1, (IQueryAtom) oldAtom2, oldBond.getOrder(), hs.getBuilder() );
			} else {
				
				try {
					newBond = (IBond) oldBond.clone();
					newBond.setAtom( oldAtom1, 0 );
					newBond.setAtom( oldAtom2, 1 );
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
			}
			
			// restore other properties
			if( oldBond.getProperty( bondFrequencyType ) != null )
				newBond.setProperty(bondFrequencyType, oldBond.getProperty( bondFrequencyType ));
			
			if( oldBond.getProperty( topologyType ) != null )
				newBond.setProperty(topologyType, oldBond.getProperty( topologyType ));
			
			if( oldBond.getProperty( bondMolOriginType ) != null )
				newBond.setProperty(bondMolOriginType, oldBond.getProperty( bondMolOriginType ));
			
			//System.err.println( "update bond " + nHs.getAtomNumber( oldAtom1 ) + " " + nHs.getAtomNumber( oldAtom2 ) );
			newBonds[b] = newBond;
			
		}
		
		
		
		
		nHs.setBonds( newBonds );

		
		return nHs;
	}


	private IQueryAtom createHSAtom( IAtom origAtom, IAtom hAtom ) {
		
		IQueryAtom newAtom = null;
		
		if( hAtom != null ) {
			// the hyperstructure atom SHOULD be a query atom
			IQueryAtom hsAtom = (IQueryAtom) hAtom;
			if( ! hsAtom.matches(origAtom) ) {
				//newAtom = new AnyAtom( DefaultChemObjectBuilder.getInstance() );
				
				SmartsAtomExpression newSMARTSAtom = null;
				if( !( hsAtom instanceof SmartsAtomExpression) ) {
					newSMARTSAtom = new SmartsAtomExpression( DefaultChemObjectBuilder.getInstance() );
					newSMARTSAtom.tokens.add( new SmartsExpressionToken(SmartsConst.AP_AtNum, PeriodicTable.getAtomicNumber( hAtom.getSymbol() ) )  );
					System.out.println( "hAtom symbol - " + hAtom.getSymbol() );
				} else {
					newSMARTSAtom = (SmartsAtomExpression) hsAtom;
				}
				
				if( ! hAtom.getSymbol().equals( origAtom.getSymbol() ) ) {
					newSMARTSAtom.tokens.add( new SmartsExpressionToken( SmartsConst.LO_OR + SmartsConst.LO, 3 )  );
					newSMARTSAtom.tokens.add( new SmartsExpressionToken(SmartsConst.AP_AtNum, PeriodicTable.getAtomicNumber( origAtom.getSymbol() ) )  );
				}
				//newSMARTSAtom.to
				
				newAtom = newSMARTSAtom;
				if( newAtom.getSymbol() == null )
					newAtom.setSymbol( "Si" );  // dummy symbol
				
			} else {
				newAtom = hsAtom;
			}
		} else if( origAtom instanceof IQueryAtom ) {
			try {
				newAtom = (IQueryAtom) origAtom.clone();
				
			} catch (CloneNotSupportedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			
		} else {
			if( ConvenienceTools.isAromatic(origAtom) ) {
				newAtom = new AromaticSymbolQueryAtom( hs.getBuilder() );
				newAtom.setSymbol( origAtom.getSymbol() );
				newAtom.setProperty( aromaticAtomProp, 5 );
				newAtom.setFlag( CDKConstants.ISAROMATIC, true );
			} else {
				newAtom = new AliphaticSymbolQueryAtom( hs.getBuilder() );
				newAtom.setSymbol( origAtom.getSymbol() );
				System.out.println( "hAtom symbol clone - " + newAtom.getSymbol() );
			}
		}
		
		// charge hack - needed for phosphorus and other atom perception 
		newAtom.setFormalCharge( 0 );
		
		//System.out.println(" new atom - " + newAtom.getSymbol() );
		
		return newAtom;
	}
	

	@Deprecated
	private IQueryBond createSMARTSRingChainBond( IAtom translated1, IAtom translated2, boolean isRing, boolean isAromatic, Order order ) {
		SmartsBondExpression sBond = new SmartsBondExpression( DefaultChemObjectBuilder.getInstance() );
		sBond.setAtom(translated1, 0);
		sBond.setAtom(translated2, 1);
		
		sBond.setFlag(CDKConstants.ISINRING, true);
		
		if( ! isRing ) {
			sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_NOT );
			sBond.setFlag(CDKConstants.ISINRING, false);
		}
		
		sBond.tokens.add( SmartsConst.BT_RING );
		
		if( order != null ) {
			sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_ANDLO );
			
			if( order == Order.SINGLE )
				sBond.tokens.add( SmartsConst.BT_SINGLE );
			
			if( order == Order.DOUBLE )
				sBond.tokens.add( SmartsConst.BT_DOUBLE );
			
			if( order == Order.TRIPLE )
				sBond.tokens.add( SmartsConst.BT_TRIPLE );
			
			sBond.setOrder( order );
		}
		
		if( isAromatic ) {
			if( order == null )
				sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_OR );
			else
				sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_ANDLO );
				
			sBond.tokens.add( SmartsConst.BT_AROMATIC );
		}
		
		
		return sBond;
	}
	
	
	/**
	 * If both bonds are in a ring (SMARTS and normal definitions), then the bond is made ring-only.
	 * If one is in a ring and the other isn't, the bond is made standard.
	 * If both bonds are not in a ring (chain bonds), then it's made chain only.
	 * 
	 * Bond order is still kept into account in this.  So is aromaticity, at least as far as rings are concerned (not chains)
	 * 
	 * @param hsBond
	 * @param origBond
	 * @param aromaticRing
	 * @return
	 */
	private IQueryBond createSMARTSRingChainBond( IAtom translated1, IAtom translated2, IBond hsBond, IBond origBond ) {
		IQueryBond newBond = null;
		Boolean bothRings = null;  // true = both in rings, false = both in chains, null = exclusive-OR
		
		
		
		if( hsBond != null && !( hsBond instanceof SmartsBondExpression) )  // if something's been marked as non-topology specific, ignore that bond
			bothRings = null;
		
		
		if( hsBond != null && ConvenienceTools.isRingBond( origBond ) && ConvenienceTools.isRingBond( hsBond ) ) { // both are in rings
			bothRings = true;
		} else if( hsBond != null && ! ConvenienceTools.isRingBond( origBond ) && ! ConvenienceTools.isRingBond( hsBond ) ) {  
			bothRings = false;
		} else if( hsBond == null  ) {  // assume that origBond is never null
			bothRings = ConvenienceTools.isRingBond( origBond );
		}
		
		//boolean isAromatic = ConvenienceTools.isAromatic(origBond) || ConvenienceTools.isAromatic(translated1) || ConvenienceTools.isAromatic(translated2);
		
		Boolean bothAromatic = null;  // true = both aromatic, false = both non-aromatic, null = make new bond map both single & aromatic
		if( hsBond != null && ConvenienceTools.isAromatic( origBond ) && ConvenienceTools.isAromatic( hsBond ) ) {
			bothAromatic = true;
		} else if( hsBond != null && ! ConvenienceTools.isAromatic( origBond ) && ! ConvenienceTools.isAromatic( hsBond ) ) {
			bothAromatic = false;
		} else if( hsBond == null && ConvenienceTools.isAromatic( origBond ) ) {  // non-existent hsBond - default to properties of other bond
			bothAromatic = true;
		} else if( ConvenienceTools.isAromatic( translated1 ) && ConvenienceTools.isAromatic( translated2 ) && ConvenienceTools.isAromatic( origBond ) ) {
			bothAromatic = true;
		} else if( ! ConvenienceTools.isAromatic( translated1 ) && ! ConvenienceTools.isAromatic( translated2 ) && ! ConvenienceTools.isAromatic( origBond ) ) {
			bothAromatic = false;
		}
		
		
		Order order = origBond.getOrder();
		if( hsBond != null && hsBond.getOrder() != origBond.getOrder() )
			order = null;
		
		/*if( bothAromatic != null && bothAromatic == true ) {
			newBond = new AromaticQueryBond( translated1.getBuilder() );
			newBond.setAtom(translated1, 0);
			newBond.setAtom(translated2, 1);
			newBond.setOrder( Order.SINGLE );
		} else*/ if( bothRings == null ) {
			
			if( bothAromatic == null ) {
				newBond = new SingleOrAromaticBond( translated1.getBuilder() );
				newBond.setAtom(translated1, 0);
				newBond.setAtom(translated2, 1);
				newBond.setOrder( Order.SINGLE );
			} else if( bothAromatic == true ) {
				newBond = new AromaticQueryBond( translated1.getBuilder() );
				newBond.setAtom(translated1, 0);
				newBond.setAtom(translated2, 1);
				newBond.setOrder( Order.SINGLE );
			} else {
				if( order != null )
					newBond = new OrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, order, hs.getBuilder() );
				else 
					newBond = new AnyOrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, Order.SINGLE, hs.getBuilder() );
			}
			
			//newBond.setFlag(CDKConstants.ISINRING, true);
			
		} else {
			System.err.println("RING1 " + origBond + "  " + hsBond );
			
			SmartsBondExpression sBond = new SmartsBondExpression( DefaultChemObjectBuilder.getInstance() );
			sBond.setAtom(translated1, 0);
			sBond.setAtom(translated2, 1);
			
			sBond.setFlag(CDKConstants.ISINRING, true);
			
			if( bothRings == false ) {
				sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_NOT );
				sBond.setFlag(CDKConstants.ISINRING, false);
				//sBond.setProperty(topologyType, "c");
			} else if (bothRings == true) {
				//sBond.tokens.add( SmartsConst.BT_AROMATIC );
				//sBond.setFlag(CDKConstants.ISAROMATIC, true);
				//sBond.setProperty(topologyType, "r");
			}
			
			sBond.tokens.add( SmartsConst.BT_RING );
				
			
			if( bothAromatic != null && bothAromatic == true ) {
				sBond.tokens.add( SmartsConst.BT_AROMATIC );
				sBond.setFlag(CDKConstants.ISAROMATIC, true);
				sBond.setOrder( Order.SINGLE );
			} else if( order != null ) {
				sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_ANDLO );
				
				if( order == Order.SINGLE )
					sBond.tokens.add( SmartsConst.BT_SINGLE );
				
				if( order == Order.DOUBLE )
					sBond.tokens.add( SmartsConst.BT_DOUBLE );
				
				if( order == Order.TRIPLE )
					sBond.tokens.add( SmartsConst.BT_TRIPLE );
				
				sBond.setOrder( order );
			}
			
			if( bothAromatic == null ) {
				//if( order == null )
					sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_OR );
				//else
					//sBond.tokens.add( SmartsConst.LO + SmartsConst.LO_ANDLO );
					
				sBond.tokens.add( SmartsConst.BT_AROMATIC );
				sBond.setFlag(CDKConstants.ISAROMATIC, true);
			}
			
			newBond = sBond;
			
		} 
		
		
		return newBond;
	}
	
	/**
	 * Creates a new bond in the hyperstructure.  Generally used when a bond doesn't exist in the MCS
	 * 
	 * @param translated1
	 * @param translated2
	 * @param origBond
	 * @param hsBond
	 * @return
	 */
	private IQueryBond createHSBond( IAtom translated1, IAtom translated2, IBond origBond, IBond hsBond ) {
		
		IQueryBond hsqBond = (IQueryBond) hsBond;  // likely to be null
		
		
		
		IQueryBond newBond = null;
		
		boolean matches = true;
		String bondTopologyType = defineBondTopology(hsBond, origBond);
		
		if( hsBond != null ) {
			if( hsqBond.getProperty( origBondType ) != null )
				hsBond = hsqBond.getProperty( origBondType );
			
			if( hsBond instanceof IQueryBond ) {
				hsqBond = (IQueryBond) hsBond;
	
				matches = hsqBond.matches(origBond);
				System.err.println("BLARG " + hsBond.getAtom(0).getSymbol() + hsBond.getAtom(1).getSymbol() + hsBond.getOrder() + " " + origBond.getOrder() + matches + "  " + hsBond + " " + origBond );
			} else {
				matches = (hsBond.getOrder() == origBond.getOrder())  ;
				System.err.println("HECK " + matches + "  " + hsBond );
			}
		}
		
		
		//System.err.println("BLA2RG  " + matches + " " + hsBond );
		
		if( ringEnforcement  ) {
			
			newBond = createSMARTSRingChainBond( translated1, translated2, hsBond, origBond );
			
		} else {
			if( ConvenienceTools.isAromatic( translated1 ) && ConvenienceTools.isAromatic( translated2 ) && ConvenienceTools.isAromatic( origBond ) ) {
				newBond = new AromaticQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), hs.getBuilder() );
			
	
			
			/*} else if( ConvenienceTools.isAromatic( translated1 ) || ConvenienceTools.isAromatic( translated2 ) ) {
				//newBond = new AromaticOrSingleQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), hs.getBuilder() );
				newBond = new SingleOrAromaticBond(hs.getBuilder() );
				newBond.setAtom(translated1, 0);
				newBond.setAtom(translated2, 1);
				newBond.setOrder(origBond.getOrder() );*/
			} else {
				if( matches ) {
					//newBond = new AnyOrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), hs.getBuilder() );
					
						if( ! (origBond instanceof IQueryBond) ) {
							
							/*if ( translated1 instanceof AtomicNumberAtom || translated2 instanceof AtomicNumberAtom ) {
								newBond = new AnyOrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), hs.getBuilder() );
							} else*/ if( ConvenienceTools.isAromatic(translated1) && ConvenienceTools.isAromatic(translated2) && ConvenienceTools.isAromatic( origBond ) ) {
								
								
								
									newBond = new AromaticQueryBond( hs.getBuilder() );
									newBond.setAtom(translated1, 0);
									newBond.setAtom(translated2, 1);
									newBond.setOrder( Order.SINGLE );
								
							// Allow aromatic bonds to map in cases where only one atom in the bond is aromatic
							} else {
								
								
									System.err.println("ORDER1 " + matches + "  " + origBond + " | " + hsBond );
									//newBond = new AnyOrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), hs.getBuilder() );
									boolean aromaticRing = ConvenienceTools.isAromatic(origBond) || ConvenienceTools.isAromatic(translated1) || ConvenienceTools.isAromatic(translated2);
									
									if( aromaticRing ) {
										newBond = new SingleOrAromaticBond( translated1.getBuilder() );
										newBond.setOrder( Order.SINGLE );
									} else {
										
										newBond = new OrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), origBond.getBuilder() );
										
									}
								
								
								newBond.setAtom(translated1, 0);
								newBond.setAtom(translated2, 1);
								
							}
							
							/*
							else if( ConvenienceTools.isAromatic(translated1) ^ ConvenienceTools.isAromatic(translated2) ) {
								
								// if both are ring-enforced bonds, make sure the new bond is also ring-enforced - @&-,:
								if( bothRings ) {
									System.err.println("RING1 " + matches + "  " + hsBond );
	
									newBond = createSMARTSRingChainBond(translated1, translated2, true, true, Order.SINGLE);
								} else if( chainOnly ) {
									System.err.println("CHAIN1 " + matches + "  " + hsBond );
									
									newBond = createSMARTSRingChainBond(translated1, translated2, false, true, Order.SINGLE);
								} else {
									newBond = new SingleOrAromaticBond( hs.getBuilder() );
								}
							
								newBond.setAtom(translated1, 0);
								newBond.setAtom(translated2, 1);
								newBond.setOrder( Order.SINGLE );
							} else {
								// if both are in rings, change to generic ring bond
								if( bothRings ) {
									System.err.println("RING2 " + matches + "  " + hsBond );
									/*
									
									newBond.setAtom(translated1, 0);
									newBond.setAtom(translated2, 1);
									newBond.setOrder( Order.SINGLE );
									*/
							/*
									newBond = createSMARTSRingChainBond(translated1, translated2, true, false, Order.SINGLE);
								} else if( chainOnly ) {
									System.err.println("CHAIN1 " + matches + "  " + hsBond );
									newBond = createSMARTSRingChainBond(translated1, translated2, false, false, Order.SINGLE);
								} else {
									newBond = new AnyOrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), hs.getBuilder() );
								}
								
							}*/
						} else {
							// clone the original Query Bond, though unlikely unless SMARTS queries are being mapped to the hyperstructure
							try {
								newBond = (IQueryBond) origBond.clone();
								newBond.setAtom(translated1, 0);
								newBond.setAtom(translated2, 1);
								//System.err.println("Blech");
							} catch (CloneNotSupportedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
				} else {
					System.err.println("MATCHES NOT");
					
					
						
						if( hsBond instanceof AromaticQueryBond && origBond.getOrder() == Order.SINGLE && ! ConvenienceTools.isAromatic( origBond ) ) {
							System.out.println("asoidjsao");
							newBond = new SingleOrAromaticBond( DefaultChemObjectBuilder.getInstance() );
							newBond.setAtom(translated1, 0);
							newBond.setAtom(translated2, 1);
							newBond.setOrder( Order.SINGLE );
						} else {
						
							newBond = new AnyOrderQueryBond( (IQueryAtom) translated1, (IQueryAtom) translated2, origBond.getOrder(), hs.getBuilder() );
						}
					
				}
			}
		}
		
		if( bondFrequencies ) {
			
			setBondFrequencyAttributes(hsBond, newBond, false);
			
			newBond.setProperty( topologyType, bondTopologyType );
			
		}
		
		
		
		return newBond;
	}
	
	
	private String defineBondTopology( IBond hsBond, IBond origBond ) {
		Boolean bothRings = null;
		
		// if the property has already been removed, don't bother setting it again
		/*if( hsBond != null && hsBond.getProperty( topologyType ) == null )
			return;*/
		
		
		/*if( 
				hsBond != null && 
				origBond.getProperty( topologyType ) != null &&
				hsBond.getProperty( topologyType ) != origBond.getProperty( topologyType ) 
		) { // remove conflicting topologies
			//hsBond.setProperty(topologyType, "p");
			return null;
		}*/
		
		
		if( hsBond != null && ConvenienceTools.isRingBond( origBond ) && ConvenienceTools.isRingBond( hsBond ) ) { // both are in rings
			bothRings = true;
		} else if( hsBond != null && ! ConvenienceTools.isRingBond( origBond ) && ! ConvenienceTools.isRingBond( hsBond ) ) {  
			bothRings = false;
		} else if( hsBond == null  ) {  // assume that origBond is never null
			bothRings = ConvenienceTools.isRingBond( origBond );
		}
		
		if( bothRings != null ) {
			if( bothRings == true ) {
				//hsBond.setProperty( topologyType, "r" );
				return "r";
			} else if( bothRings == false ) {
				//hsBond.setProperty( topologyType, "c" );
				return "c";
			} 
		} else {
			//hsBond.removeProperty( topologyType );
		}
		
		return null;

	}
	
	
	
	/**
	 * Update bond origin information (and frequency information) for a hyperstructure
	 * 
	 * @param hsBond
	 * @param newBond
	 * @param increment
	 */
	private void setBondFrequencyAttributes( IBond hsBond, IBond newBond, boolean increment ) {
		
		if( hsBond != null /*&& hsBond.getProperty(bondFrequencyType) != null*/ ) {  
			
			if( increment ) {
				int freq = (Integer) hsBond.getProperty(bondFrequencyType);
				newBond.setProperty(bondFrequencyType, freq + 1 );
				hsBond.setProperty(bondFrequencyType, freq + 1 );
			} else {
				newBond.setProperty(bondFrequencyType, (Integer) hsBond.getProperty(bondFrequencyType) );
			}
		} else {
			newBond.setProperty(bondFrequencyType, 1 );
			//System.err.println("HSB null freq");
		}
		

		//newBond.setProperty( topologyType, bondTopologyType );

		if( queryMolId > 0 ) {
			List<Integer> molOrigins = null;
			
			if( hsBond != null && hsBond.getProperty(bondMolOriginType) != null ) {
				molOrigins = (List<Integer>) hsBond.getProperty(bondMolOriginType);
				
				if( increment )
					molOrigins.add( queryMolId );
				//System.out.println("molOrigins - " + molOrigins);
			} else {
				molOrigins = new ArrayList<Integer>(10);
				molOrigins.add( queryMolId );
			}
			
			
			newBond.setProperty(bondMolOriginType, molOrigins);
		}
	}

	
	public Map<IAtom, IAtom> atomMapFromChromosome( List<Integer> match ) {
		
		HashMap<IAtom, IAtom> atomMap = new HashMap<IAtom, IAtom>();

	 	   // translate mapping from array to hash
	 	   for( int m = 0; m < match.size(); m++ ) {
	 		   if( match.get(m) >= 0 ) {
	 			   if( hs.getAtom( m ) != null && queryMol.getAtom( match.get(m) ) != null )
	 				   atomMap.put( hs.getAtom( m ), queryMol.getAtom(  match.get(m) ) );
	 		   }
	 	   }
	 	   
	 	   return atomMap;
	}
	
	
	public IQueryAtomContainer createHyperstructure( ArrayList<Integer> match ) {
	 	   
	 	   return createHyperstructure( atomMapFromChromosome(match) );
	}

	/**
	 * Finds out bonds in current mapping that have been mapped
	 * 
	 * Append any bonds that are not in the mapping to the hyperstructure.  For all bonds in the query, if the bond is not in the mapping then add it to the hyperstructure.  
	 * 
	 * To find a bond in the hyperstructure from the query, get the 2 atoms from the bond, put into the atom map, then see if there's a bond between the 2 translated atoms in the hyperstructure.
	 * 
	 * If bond not found:
	 * - both atoms exist but there's no bond - add a bond
	 * - both atoms exist but with a different bond - add degenerate bond type
	 * - one of the atoms doesn't exist - add the "null" atom to the hyperstructure and connect with a bond.  Add the new atom as a corresponding pair to the supplied "mapping" argument
	 * - both atoms don't exist - add them to the hyperstructure & connect them, adding both atoms to the mapping argument.
	 * 
	 * In the case of adding to the mapping argument, this means that any unused atoms we look through later on can be reconnected.
	 * 
	 * @param mapping
	 * @return
	 */
	public IQueryAtomContainer createHyperstructure( Map<IAtom, IAtom> mapping ) {

		
		/* 
		 * Construct a hash which translates the "old" index from the chromosome to the newly assigned
		 * index in the hyperstructure, for an unmapped atom that's been added to the hyperstructure
		 */
		
		Map<IBond, IBond> mappedBonds = ConvenienceTools.makeBondMapOfAtomMap(queryMol, hs, mapping, false); // bad practice, not the same as the actual bond mapping
		
		return createHyperstructure( mapping, mappedBonds );
		
	}
	
	
	public IQueryAtomContainer createHyperstructure( Map<IAtom, IAtom> mapping, Map<IBond, IBond> mappedBonds ) {	
		
		IQueryAtomContainer novelHs = null;
		novelHs = (IQueryAtomContainer) hs;
		ConvenienceTools.countRings( novelHs );
		

		/*if( bondFrequencies ) {
			for( IBond b : novelHs.bonds() ) {
				//b.setProperty(bondFrequencyType, 1);
				IBond orig = (IBond) b.getProperty( origBondType );
				orig.setProperty(bondFrequencyType, 1);
			}
		}*/
		
		//HashMap<IAtom, Integer> usedAtoms = new HashMap<IAtom, Integer>(); 
		
		
		// go through all bonds in the query - adding those bonds (with atoms) that don't exist in the hyperstructure
		for( IBond b : queryMol.bonds() ) {
			IAtom translated1 = null, translated2 = null; 
			if( mapping.containsKey( b.getAtom(0) ) )
				translated1 = mapping.get( b.getAtom(0) );
			 
			if( mapping.containsKey( b.getAtom(1) ) )
				translated2 = mapping.get( b.getAtom(1) );
			
			if( ! mappedBonds.containsKey(b) ) {  // bond not present in MCS
				
				//System.out.println( "translated atoms: " + translated1 + " " + translated2 );
				
				IBond hsBond = null;
				
				
				//if( translated1 != null && translated2 != null ) { // use the atom mappings to attempt to find a bond in the hyperstructure
				//	hsBond = novelHs.getBond(translated1, translated2);
				//}	
					/*
					 * Special case where one bond doesn't match another (as opposed to a series of unmatched bonds in a row)
					 * 
					 * non-ring bond
					 * Choose one atom on the bond in the hyperstructure, and an atom on an adjacent bond which passes this unmatching bond.  
					 * Create 
					 */
					/*if( hsBond != null && hsBond.getOrder() != b.getOrder() ) {
						//SMARTSBond hsSMARTSbond = (SMARTSBond) hsBond;
						
						//if( ! hsSMARTSbond.matches(b) ) {
							translated1 = createHSAtom( b.getAtom(0), null );
							translated2 = createHSAtom( b.getAtom(1), null );
							
							novelHs.addAtom( translated1 );
							novelHs.addAtom( translated2 );
							
							mapping.put( b.getAtom(0), translated1 );  // this is needed in-case the same atom is referred to in subsequent query bonds
							mapping.put( b.getAtom(1), translated2 );  // this is needed in-case the same atom is referred to in subsequent query bonds
							
							IQueryBond newBond = createHSBond( translated1, translated2, b, mappedBonds.get(b) );
							novelHs.addBond( newBond );
							
							mappedBonds.remove( b );
							mappedBonds.put( b, novelHs.getBond(translated1, translated2) );
							hsBond = null;
						//}
					}
				}
				*/
				
				if( translated1 == null ) {
					translated1 = createHSAtom( b.getAtom(0), null );
					
					mapping.remove( b.getAtom(0) );
					mapping.put( b.getAtom(0), translated1 );  // this is needed in-case the same atom is referred to in subsequent query bonds
					
					novelHs.addAtom( translated1 );
				}
					
					
				if( translated2 == null ) {
					
					translated2 = createHSAtom( b.getAtom(1), null );
				
					mapping.remove( b.getAtom(1) );
					mapping.put( b.getAtom(1), translated2 );  // this is needed in-case the same atom is referred to in subsequent query bonds
					
					novelHs.addAtom( translated2 );
				}
				
				/*
				// deal with aromatic-to-non aromatic mappings
				if( ConvenienceTools.isAromatic( b.getAtom(0) ) && ConvenienceTools.isAromatic( translated1 ) ) {
					translated1.setProperty( aromaticAtomProp, 5 );
				} else if( ConvenienceTools.isAromatic( b.getAtom(0) ) ^ ConvenienceTools.isAromatic( translated1 ) ) {
					translated1.setProperty( atomicNumberProp, 5 );
				}
				
				if( ConvenienceTools.isAromatic( b.getAtom(1) ) && ConvenienceTools.isAromatic( translated2 ) ) {
					translated2.setProperty( aromaticAtomProp, 5 );
				} else if( ConvenienceTools.isAromatic( b.getAtom(1) ) ^ ConvenienceTools.isAromatic( translated2 ) ) {
					translated2.setProperty( atomicNumberProp, 5 );
				}
				*/
				
				// case where no bond exists
				if( hsBond == null ) {
					if( translated1 != null && translated2 != null ) { // use the atom mappings to attempt to find a bond in the hyperstructure
						hsBond = novelHs.getBond(translated1, translated2);
					}	
					/*
					translated1 = createHSAtom( b.getAtom(0), null );
					translated2 = createHSAtom( b.getAtom(1), null );
					
					novelHs.addAtom( translated1 );
					novelHs.addAtom( translated2 );
					
					mapping.put( b.getAtom(0), translated1 );  // this is needed in-case the same atom is referred to in subsequent query bonds
					mapping.put( b.getAtom(1), translated2 );  // this is needed in-case the same atom is referred to in subsequent query bonds
					*/
					IQueryBond newBond = createHSBond( translated1, translated2, b, hsBond );
					novelHs.addBond( newBond );
					System.err.println( "newbond - " + newBond );
					mappedBonds.remove( b );
					mappedBonds.put( b, newBond );
					
					//defineBondTopology( newBond, b );
					System.out.println( "new bond - " + novelHs.getAtomNumber( newBond.getAtom(0) ) + " " + novelHs.getAtomNumber( newBond.getAtom(1) ) );
					
				}
				
				
				
				// a bond does exist
			} else { //else if( translated1 != null && translated2 != null && hsBond != null ) { // a bond exists
				//novelHs.getBond(translated1, translated2).setOrder( Order.QUADRUPLE );
				
				
				IBond hsBond = mappedBonds.get(b);
				IQueryBond hsqBond = (IQueryBond) hsBond;
				
				
				
				if( hsqBond.getProperty( origBondType ) != null )
					hsBond = hsqBond.getProperty( origBondType );
				
				String bondTypeTopology = defineBondTopology( hsBond, b );
				
				
				IQueryBond newBond = null;
				IQueryAtom tr1 = (IQueryAtom) hsBond.getAtom(0);
				IQueryAtom tr2 = (IQueryAtom) hsBond.getAtom(1);
				System.err.println("BAD2 " + "  " + hsBond + " " + hsBond.getProperty( origBondType ) + " " + tr1 + " " + tr2 );
				
				boolean matches = true;
				
				if( hsBond instanceof IQueryBond ) {
					hsqBond = (IQueryBond) hsBond;
					//System.out.println( "type: " + hsqBond );
					
					matches = hsqBond.matches(b);
					System.err.println("BAD " + matches + "  " + b + " " + hsqBond );
				
				} else {
					matches = hsBond.getOrder() == b.getOrder()  ;
					System.err.println("NORM " + matches + "  " + hsBond );
				}
				
				
				if( ConvenienceTools.isAromatic( b.getAtom(0) ) && ConvenienceTools.isAromatic( tr1 ) ) {
					tr1.setProperty( aromaticAtomProp, 5 );
				} else if( ConvenienceTools.isAromatic( b.getAtom(0) ) ^ ConvenienceTools.isAromatic( tr1 ) ) {
					tr1.setProperty( atomicNumberProp, 5 );
				}
				
				if( ConvenienceTools.isAromatic( b.getAtom(1) ) && ConvenienceTools.isAromatic( tr2 ) ) {
					tr2.setProperty( aromaticAtomProp, 5 );
				} else if( ConvenienceTools.isAromatic( b.getAtom(1) ) ^ ConvenienceTools.isAromatic( tr2 ) ) {
					tr2.setProperty( atomicNumberProp, 5 );
				}
				
				
				if( bondFrequencies ) {
					// update frequency
					//hsBond.setProperty(bondFrequencyType, (Integer) hsBond.getProperty(bondFrequencyType) + 1 );
					//setBondFrequencyAttributes(hsBond, hsBond, true);
					setBondFrequencyAttributes(novelHs.getBond( tr1, tr2 ), novelHs.getBond( tr1, tr2 ), true);
					
					System.err.println( "frec: " +  hsBond.getProperty(bondFrequencyType) + " , " +  hsBond.getProperty(topologyType) );
					//defineBondTopology(hsBond, b);
				}
				
				
				if( hsBond instanceof AnyOrderQueryBond ) {
					
					hsBond.setProperty(manyOrderProp, 5);

					
					
					// bonds don't match - so make a SMARTS-like bond to account for both bonds
				} else if( ! matches ) {
					
					//novelHs.removeBond( mappedBonds.get(b) );  // be careful - doens't complain if the bond does not exist!
					novelHs.removeBond( hsBond );

					newBond = createHSBond(tr1, tr2, b, hsBond);
						
						
						
					novelHs.addBond( newBond );
					
					mappedBonds.remove( b );
					mappedBonds.put( b, newBond );
					
					System.err.println( "added: " + b.getAtom(0).getSymbol() + b.getAtom(1).getSymbol() + hsBond.getOrder() + "  " + newBond );
				
				} else if( hsBond.getClass() == b.getClass() ||  matches ) {  // keep the bond (not really needed except for debugging)
					
					System.err.println( "kept: " + b.getAtom(0).getSymbol() + b.getAtom(1).getSymbol() + hsBond.getOrder() + "  " + hsBond + " | " + b );
					/*
					novelHs.removeBond( hsBond );
					//new SingleOrAromaticBond(hs.getBuilder()).;
					newBond = new AromaticOrSingleQueryBond( tr1, tr2, b.getOrder(), hs.getBuilder() );
					novelHs.addBond( newBond );*/
					 
				} 
				
				
				
				// Deal with ring inconsistencies - notably null bonds
				if( ringEnforcement && hsBond instanceof SmartsBondExpression ) {
					//ConvenienceTools.countRings( novelHs );
					IBond hsBondR = novelHs.getBond( tr1, tr2 );
					
					// hack, sometimes hsBondR is null
					if( hsBondR == null )
						hsBondR = novelHs.getBond( translated1, translated2 );
					
					//boolean aromaticRing = ConvenienceTools.isAromatic(b) || ConvenienceTools.isAromatic(tr1) || ConvenienceTools.isAromatic(tr2);
					
					
					// hack, sometimes hsBondR is null
					novelHs.removeBond( hsBondR );
					
					
					newBond = createSMARTSRingChainBond( translated1, translated2, hsBond, b );
					//defineBondTopology( newBond, b );
					//System.err.println( "new ring bond: " + newBond );
					
					if( bondFrequencies ) {
						// update frequency
						//newBond.setProperty(bondFrequencyType, hsBond.getProperty(bondFrequencyType) );
						setBondFrequencyAttributes(hsBond, newBond, false);
						//newBond.setProperty(topologyType, bondTypeTopology);
						//defineBondTopology( hsBond, b );
					}
					
					System.err.println( "rings: " + tr1.getSymbol() + tr2.getSymbol() + hsBond.getOrder() + " " + b.getAtom(0).getSymbol() + b.getOrder() + "  " + ConvenienceTools.isRingBond( hsBond ) + " | " + ConvenienceTools.isRingBond( b ) + hsBond + " , " + newBond.getProperty(topologyType) );
					
					
					novelHs.addBond( newBond );
					
					mappedBonds.remove( b );
					mappedBonds.put( b, newBond );
				}
				
				// set topology type
				//System.err.println( novelHs.getBond( tr1, tr2 ) );
				if( novelHs.getBond( tr1, tr2 ) == null ) {
					boolean lolss = true;
				}
				
				novelHs.getBond( tr1, tr2 ).setProperty(topologyType, bondTypeTopology);
				//hsBond.setProperty(topologyType, bondTypeTopology);
				
				//System.err.println( "frec2: " +  novelHs.getBond( tr1, tr2 ).getProperty(bondFrequencyType) + " , " +  novelHs.getBond( tr1, tr2 ).getProperty(topologyType) );
				
					//novelHs.addAtom( tr1 );
					//novelHs.addAtom( tr2 );
				
				/*
				mapping.remove( b.getAtom(0) );
				mapping.remove( b.getAtom(1) );
				mapping.put( b.getAtom(0), tr1 );
				mapping.put( b.getAtom(1), tr2 );
				*/
				//mappedBonds.put( b, novelHs.getBond(tr1, tr2) );
				
				
			}
			/*
			IAtom[] newAtoms = new IAtom[ novelHs.getAtomCount() ] ;
			HashMap<IAtom, IAtom> oldToNewAtoms = new HashMap<IAtom, IAtom>();
			
			for( int a = 0; a < novelHs.getAtomCount(); a++ ) {
				newAtoms[a] = novelHs.getAtom(a);
			}
			
			for( int a = 0; a < queryMol.getAtomCount(); a++ ) {
				
				IQueryAtom hsAtom = (IQueryAtom) mapping.get( queryMol.getAtom(a) );
				
				if( hsAtom != null ) {
					if( ! hsAtom.matches( queryMol.getAtom(a) ) ) {
						int index = novelHs.getAtomNumber( hsAtom );
						if( index >= 0 ) {
							IAtom newAtom = createHSAtom( queryMol.getAtom(a), hsAtom );
							System.out.println( "the new atom is " + newAtom + " " + index + " " + a );
							
							oldToNewAtoms.put(hsAtom, newAtom);
							newAtoms[index] = newAtom;
							mapping.put( queryMol.getAtom(a), newAtom );
						}
					}
				}
			}
			
			for( IBond hb : novelHs.bonds() ) {
				if( oldToNewAtoms.containsKey( hb.getAtom(0) )  ) {
					hb.setAtom( oldToNewAtoms.get(hb.getAtom(0)), 0);
					
				}
				if( oldToNewAtoms.containsKey( hb.getAtom(1) )  )
					hb.setAtom( oldToNewAtoms.get(hb.getAtom(1)), 1);
				
				if( hb.getProperty(origBondType) != null ) {
					IBond origBond = hb.getProperty(origBondType);
					
					if( oldToNewAtoms.containsKey( origBond.getAtom(0) )  ) {
						origBond.setAtom( oldToNewAtoms.get(origBond.getAtom(0)), 0);
						System.out.println( "bond replaced " );
					}
					
					if( oldToNewAtoms.containsKey( origBond.getAtom(1) )  )
						origBond.setAtom( oldToNewAtoms.get(origBond.getAtom(1)), 1);
				}
			}
			
			novelHs.setAtoms(newAtoms);
			*/
			
			/*
			// go through all bonds in the query - adding those bonds (with atoms) that don't exist in the hyperstructure
			for( Entry<IAtom, IAtom> e : mapping.entrySet() ) {
				if( e.getKey() != null && e.getValue() != null ) {
					IQueryAtom hsAtom = (IQueryAtom) e.getValue();
					
					if( ! hsAtom.matches( e.getKey() ) ) {
						IAtom newAtom = createHSAtom(e.getKey(), hsAtom);
						newAtom.setProperty("new", "new");
						int index = novelHs.getAtomNumber( e.getValue() );
						if( index >= 0 ) {

							for( IBond cb : novelHs.getConnectedBondsList( hsAtom ) ) {
								if( cb.getAtom(0) == hsAtom  )
									cb.setAtom(newAtom, 0);
								
								if( cb.getAtom(1) == hsAtom  )
									cb.setAtom(newAtom, 1);
								
								System.out.println( "the new atom is " + newAtom + " " + index );
							}
							
							
							novelHs.setAtom(index, newAtom);
							
						}
					}
				}
				
				
			}*/
			//IBond hsBond = mappedBonds.get(b);
			//System.err.println( "frec: " +  hsBond.getProperty(bondFrequencyType) + " , " +  hsBond.getProperty(topologyType) );
		}
		
		/*if( bondFrequencies ) {
			for( IBond hb : novelHs.bonds() ) {
				if( hb instanceof SmartsBondExpression ) {
					SmartsBondExpression sBond = (SmartsBondExpression) hb;
					
					if(    sBond.tokens.get(0) == SmartsConst.LO + SmartsConst.LO_ANDLO 
						&& sBond.tokens.get(1) == SmartsConst.BT_RING ) {
						sBond.tokens.remove(1);
						sBond.tokens.remove(0);
						System.out.println("locasdhut");
					} else if( sBond.tokens.get(0) == SmartsConst.BT_RING ) {
						sBond.tokens.remove(0);
					}
				}
			}
		}*/
		
		
		
		//ConvenienceTools.countRings( novelHs );
		
		return novelHs;
	}
	
	public double fitness(ArrayList<Integer> mapChr) {
		
		HashMap<IAtom, IAtom> atomMap = new HashMap<IAtom, IAtom>();
	 	   
	 	   //System.out.println("stuff: " + mapChr.size() + " " + mapper.getQueryMol().getAtomCount() + " " + mapper.getMainMol().getAtomCount() + " | " + mapChr );
	 	   
	 	   // translate mapping from array to hash
	 	   for( int m = 0; m < mapChr.size(); m++ ) {
	 		   if( mapChr.get(m) >= 0 ) {
	 			   if( hs.getAtom( mapChr.get(m) ) != null && queryMol.getAtom( m ) != null )
	 				   atomMap.put( hs.getAtom( mapChr.get(m) ), queryMol.getAtom(  m ) );
	 		   }
	 	   }
			
	 	  Map<IBond, IBond> bondMap = ConvenienceTools.makeBondMapOfAtomMap( hs, queryMol, atomMap );
	 	 //System.out.println( bondMap.size() + " " + bondMap );
	 	  return bondMap.size() / (double) queryMol.getBondCount() ;
		}
	
	//protected IQueryAtomContainer hs;
	//protected IAtomContainer queryMol;
	private String manyOrderProp = "_order";
	private String atomicNumberProp = "_an";
	private String aromaticAtomProp = "_ar";
	private String origBondType = "_obt";
	
	public static final String bondFrequencyType = "_bfreq";
	public static final String bondMolOriginType = "_bMOrigin";
	public static final String topologyType = "_topType";
	
	private boolean ringEnforcement = false;
	private boolean bondFrequencies = true;
	private int queryMolId;
}
