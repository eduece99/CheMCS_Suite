package org.cisrg.knime;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.cisrg.ambit.SmartsHelper;
import org.cisrg.ambit.SmartsParser;
import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.MCSMethods;
import org.knime.base.node.parallel.appender.AppendColumn;
import org.knime.base.node.parallel.appender.ColumnDestination;
import org.knime.base.node.parallel.appender.ExtendedCellFactory;
import org.knime.base.node.parallel.appender.ThreadedColAppenderNodeModel;
import org.knime.chem.types.SmartsCell;
import org.knime.chem.types.SmartsValue;
import org.cisrg.hyperstructures.CDKSMARTSHyperstructureFitness;
import org.knime.core.data.AdapterValue;
import org.knime.core.data.DataCell;
import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataRow;
import org.knime.core.data.DataTable;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataType;
import org.knime.core.data.collection.CollectionCellFactory;
import org.knime.core.data.collection.ListCell;
import org.knime.core.data.def.DoubleCell;
import org.knime.core.data.def.IntCell;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.knime.commons.CDKNodeUtils;
import org.openscience.cdk.knime.type.CDKCell3;
import org.openscience.cdk.knime.type.CDKValue;
import org.openscience.cdk.renderer.generators.standard.StandardGenerator;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;


/**
 * This is the model implementation of GraphSimilarity.
 * Graph-based similarity searching with support for hyperstructures
 *
 * @author CISRG (Edmund Duesbury)
 */
public class GraphSimilarityNodeModel extends ThreadedColAppenderNodeModel {
    
    // the logger instance
    private static final NodeLogger logger = NodeLogger
            .getLogger(GraphSimilarityNodeModel.class);
        
    /** the settings key which is used to retrieve and 
        store the settings (from the dialog or from a settings file)    
       (package visibility to be usable from the dialog). */
    /** the settings key which is used to retrieve and 
    store the settings (from the dialog or from a settings file)    
   (package visibility to be usable from the dialog). */
	static final String CFGKEY_REFMOLCOLUMN = "_refMolColumn";
	static final String CFGKEY_DATABASEMOLCOLUMN = "_databaseMolColumn";
	static final String CFGKEY_ALGORITHM = "_mappingAlgorithm";
	static final String CFGKEY_AGGREGATION = "_aggregation";
	static final String CFGKEY_RETURNTYPE = "_returnType";
	static final String CFGKEY_SIMILARITYTYPE = "_similarityType";
	static final String CFGKEY_OUTPUTMCS = "_outputMCS";
	static final String CFGKEY_MOREGHOSTINFO = "_moreGhostInfo";
	static final String CFGKEY_TOPODISTANCECONSTRAINT = "_topoDistanceConstraint";
	static final String CFGKEY_RAYMONDHEURISTICS = "_raymondHeuristics";
	static final String CFGKEY_RINGHEURISTICS = "_ringHeuristics";
	static final String CFGKEY_EXPTIMEOUT = "_expansionTimeOut";
	
	
	static final String similarityColumnName = "Similarity";
	static final String tanimotoColumnName = "Tanimoto"; 
	static final String tverskyColumnName = "Tversky";
	static final String commonBondNumberColumnName = "Common Bonds";
	static final String commonBondWeightsColumnName = "Common Bond Weights";
	static final String commonBondTopologiesColumnName = "Common Bond Topologies";
	static final String refBondNumberColumnName = "Reference Bonds";
	static final String refWeightedSumColumnName = "Reference Weighted Sum";
	static final String databaseBondNumberColumnName = "Database Bonds";
	static final String mcsFragmentSizesColumnName = "MCS fragment sizes";
	static final String mcsTimeName = "MCS Search Time";
	static final String mcsModProdTimeName = "ModProd Construction Time";
	static final String mcsModProdNodeCountName = "ModProd Nodes";
	static final String mcsModProdDensityName = "ModProd Density";
	static final String mcsUniqueSubstructures = "MCS unique substructures";
	static final String mcsUniqueGhostSubstructures = "MCS unique ghost substructures";
	static final String mcsUniqueGhostSubstructuresBondIdsName = "MCS ghost substructure Bond Indices";
	static final String mcsUniqueGhostSubstructuresSMARTSName = "MCS ghost substructure SMARTS";
	static final String referenceColumnName = "Reference";
	static final String mcsColumnName = "MCS";
    
/*	public enum MappingAlgorithmName {
		MCSPlus,
		SubStructure,
		CDKMCS,
		VFLib_cMCS,
		//vfLibGA_dMCS,
		Zhu_AER_dMCS,
		ChemAxon_cMCS,
		ChemAxon_dMCS
	};*/
	
	public static ExtendedAlgorithm[] MappingAlgorithmNames = new ExtendedAlgorithm[] {
		ExtendedAlgorithm.DEFAULT,
		ExtendedAlgorithm.CDKMCS,
		ExtendedAlgorithm.VFLibMCS,
		ExtendedAlgorithm.MCSPlus,
		ExtendedAlgorithm.ChemAxon_cMCES,
		ExtendedAlgorithm.ChemAxon_dMCES,
		ExtendedAlgorithm.consR_dMCES,
		ExtendedAlgorithm.BK_dMCES,
		ExtendedAlgorithm.BK_cMCES,
		ExtendedAlgorithm.CP_dMCES,
		ExtendedAlgorithm.RASCAL_dMCES,
		ExtendedAlgorithm.Depolli_dMCES,
		ExtendedAlgorithm.kCombu_dMCES,
		ExtendedAlgorithm.kCombu_cMCES,
		ExtendedAlgorithm.fMCS
	};
	
	
	/** Enum for the different aggregation methods. */
	public enum AggregationMethod {
		Minimum, Maximum, Average, Matrix
	}
	
	/** Enum for the different return types. */
	public enum ReturnType { 
		String, Collection
	}
	
	/** Enum for the similarity measure types. */
	public enum SimilarityType { 
		Tanimoto, Tversky, Number
	}

	
	private ExtendedAlgorithm mappingAlgorithmName = MappingAlgorithmNames[1];
	private AggregationMethod aggregationMethod = AggregationMethod.Maximum;
	private ReturnType returnType = ReturnType.String;  // default value to shut up initial configuration stages of node
	private SimilarityType simType = SimilarityType.Tanimoto;
	private String refMolColumn;
	private String databaseMolColumn;
    private boolean outputMCS;
    private boolean useWeights;
    private boolean ghostSubstructureInfo;
    private boolean detailedGhostInfo = false;
    private int topologicalDistanceLimit = -1;
    private boolean useRaymondHeuristics = false;
    private boolean useRingHeuristics = false;
    private int expansionTimeOut = 5000;
    
    private int rowCount;
    private long startTime = 0;
    private int tverskyCIndex = 0, 
    		tanimotoCIndex = 1,
    		commonBondsCIndex = 2,
    		refBondsCIndex = 3,
    		refWeightSumCIndex = 4,
    		dbBondsCIndex = 5,
    		mcsWeightsCIndex = 6,
    		mcsTopologiesCIndex = 7,
    		fragmentSizesCIndex = 8,
    		mcsTimeCIndex = 9,
    		modProdTimeCIndex = 10,
    		modProdSizeCIndex = 11,
    		modProdDensityCIndex = 12,
    		uniqueSubstructuresCIndex = 13,
    		uniqueGhostSubstructuresCIndex = 14,
    		refNamesCIndex = 15,
    		outputMCSIndex = 16, 
    		ghostSMARTSIndex = 17, 
    		ghostBondIDsIndex = 18;
    //private ExtendedIsomorphism similarityComparator = null;
    //private MCSMethods mapper = null;
    //private SmartsParser sParser = null;
	private SmartsHelper smaH;
	private SmilesGenerator sGenerator;
    
    
	/* 
	if( commonBonds > 0 ) {
		try {
			cells[6] = CollectionCellFactory.createListCell( mcsWeights );
			cells[7] = CollectionCellFactory.createListCell( mcsTopologies );
			
			cells[9] = new IntCell( modProdTime );
			cells[10] = new IntCell( modProdSize );
			cells[11] = new DoubleCell( modProdDensity );
			
			cells[12] = new IntCell( uniqueSubstructures );
			cells[13] = new IntCell( uniqueGhostSubstructures );*/
    
    /**
     * Constructor for the node model.
     */
    protected GraphSimilarityNodeModel() {
    
        // one incoming port and one outgoing port is assumed
        super(2, 1);
        
        //sParser = new SmartsParser();
        smaH = new SmartsHelper( DefaultChemObjectBuilder.getInstance() );
        sGenerator = new SmilesGenerator();
        
        int numThreads = Math.max(1, CDKNodeUtils.getMaxNumOfThreads() - 1 );  // want at least 1 core remaining to not hijack the GUI!
        this.setMaxThreads( numThreads );
        //this.setMaxThreads( 7 );
    }
    
    
    private DataColumnSpec[] createSpec(final DataTableSpec oldSpec) {

		DataColumnSpec[] outSpec = null;
		ArrayList<DataColumnSpec> tempSpec = new ArrayList<DataColumnSpec>();
		
		/*
		if ( aggregationMethod == AggregationMethod.Average) {
			DataColumnSpec colSpec = new DataColumnSpecCreator( similarityColumnName, DoubleCell.TYPE).createSpec();
			outSpec = new DataColumnSpec[] { colSpec };
		} else if ( aggregationMethod == AggregationMethod.Matrix) {
			DataColumnSpec colSpec = new DataColumnSpecCreator( similarityColumnName, ListCell.getCollectionType(DoubleCell.TYPE))
					.createSpec();
			outSpec = new DataColumnSpec[] { colSpec };
		} else {
			DataColumnSpec colSpec1 = new DataColumnSpecCreator( similarityColumnName, DoubleCell.TYPE).createSpec();
			DataColumnSpec colSpec2 = null;
			if ( returnType.equals(ReturnType.String)) {
				colSpec2 = new DataColumnSpecCreator( referenceColumnName, StringCell.TYPE).createSpec();
			} else if (returnType.equals(ReturnType.Collection)) {
				colSpec2 = new DataColumnSpecCreator( referenceColumnName, ListCell.getCollectionType(StringCell.TYPE))
						.createSpec();
			}
			outSpec = new DataColumnSpec[] { colSpec1, colSpec2 };
		}
		 */
		
		
		ArrayList<DataColumnSpec> colSpecs = new ArrayList<DataColumnSpec>(9);
		if ( aggregationMethod == AggregationMethod.Matrix) {
			colSpecs.add( new DataColumnSpecCreator( similarityColumnName, ListCell.getCollectionType(DoubleCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( tanimotoColumnName, ListCell.getCollectionType(DoubleCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( commonBondNumberColumnName, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( refBondNumberColumnName, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( refWeightedSumColumnName, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( databaseBondNumberColumnName, IntCell.TYPE ).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( commonBondWeightsColumnName, ListCell.getCollectionType(
					ListCell.getCollectionType( IntCell.TYPE) 
			) ).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( commonBondTopologiesColumnName, ListCell.getCollectionType(
					ListCell.getCollectionType( StringCell.TYPE) 
			) ).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsFragmentSizesColumnName, ListCell.getCollectionType(
					ListCell.getCollectionType( IntCell.TYPE ) 
			) ).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsTimeName, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsModProdTimeName, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsModProdNodeCountName, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsModProdDensityName, ListCell.getCollectionType(DoubleCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsUniqueSubstructures, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsUniqueGhostSubstructures, ListCell.getCollectionType(IntCell.TYPE)).createSpec() );
		} else {
			colSpecs.add( new DataColumnSpecCreator( similarityColumnName, DoubleCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( tanimotoColumnName, DoubleCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( commonBondNumberColumnName, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( refBondNumberColumnName, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( refWeightedSumColumnName, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( databaseBondNumberColumnName, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( commonBondWeightsColumnName, ListCell.getCollectionType(
					IntCell.TYPE 
			) ).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( commonBondTopologiesColumnName, ListCell.getCollectionType(
					StringCell.TYPE 
			) ).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsFragmentSizesColumnName, ListCell.getCollectionType(
				 IntCell.TYPE 
			) ).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsTimeName, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsModProdTimeName, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsModProdNodeCountName, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsModProdDensityName, DoubleCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsUniqueSubstructures, IntCell.TYPE).createSpec() );
			colSpecs.add( new DataColumnSpecCreator( mcsUniqueGhostSubstructures, IntCell.TYPE).createSpec() );
		}
		
		for( DataColumnSpec colSpec : colSpecs )
			tempSpec.add(colSpec);
		
		
		DataColumnSpec colSpec4 = null;
		if ( returnType.equals(ReturnType.String)) {
			colSpec4 = new DataColumnSpecCreator( referenceColumnName, StringCell.TYPE).createSpec();
		} else if (returnType.equals(ReturnType.Collection)) {
			colSpec4 = new DataColumnSpecCreator( referenceColumnName, ListCell.getCollectionType(StringCell.TYPE))
					.createSpec();
		}
		
		if( colSpec4 != null )
			tempSpec.add(colSpec4);
		
		if( outputMCS ) {
				if ( aggregationMethod == AggregationMethod.Matrix) {
					tempSpec.add( new DataColumnSpecCreator( mcsColumnName, ListCell.getCollectionType(
							SmartsCell.TYPE 
						) ).createSpec() );
				} else {
					tempSpec.add( new DataColumnSpecCreator( mcsColumnName, CDKCell3.TYPE ).createSpec() );
				}
		}
		
		if( detailedGhostInfo ) {
			if ( aggregationMethod == AggregationMethod.Matrix) {
				/*tempSpec.add( new DataColumnSpecCreator( mcsUniqueGhostSubstructuresBondIdsName, ListCell.getCollectionType(
						SmartsCell.TYPE 
					) ).createSpec() );*/
				tempSpec.add( new DataColumnSpecCreator( mcsUniqueGhostSubstructuresSMARTSName, ListCell.getCollectionType(
						 SmartsCell.TYPE 
					) ).createSpec() );
			} else {
				/*tempSpec.add( new DataColumnSpecCreator( mcsUniqueGhostSubstructuresBondIdsName, ListCell.getCollectionType(
						IntCell.TYPE 
				) ).createSpec() );*/
				tempSpec.add( new DataColumnSpecCreator( mcsUniqueGhostSubstructuresSMARTSName,  SmartsCell.TYPE ).createSpec() );
				
			}
		}
		
		//outSpec = new DataColumnSpec[] { colSpec1, colSpec2 };
		outSpec = tempSpec.toArray( new DataColumnSpec[0] );
		
		return outSpec;
	}



    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO Code executed on reset.
        // Models build during execute are cleared here.
        // Also data handled in load/saveInternals will be erased here.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {
        
        // check if user settings are available, fit to the incoming
        // table structure, and the incoming types are feasible for the node
        // to execute. If the node can execute in its current state return
        // the spec of its output data table(s) (if you can, otherwise an array
        // with null elements), or throw an exception with a useful user message
    	
    	/*
    	 * Auto-set the database and reference molecule settings if manual configuration has failed for some reason.
    	 * Do this by finding compatible columns for each input port.
    	 */
    	
    	System.err.println( "Configuration started "  );
    	logger.warn( "Configuration started "  );
    	
    	// (re)set "hidden" variables
    	ghostSubstructureInfo = false;
    	useWeights = false;
    	
    	if ( databaseMolColumn == null
				|| (inSpecs[0].findColumnIndex(databaseMolColumn)) == -1) {
			String dName = null;
			for (DataColumnSpec s : inSpecs[0]) {
				if (s.getType().isCompatible(CDKValue.class)) {
					dName = s.getName();
				}
			}
			if (dName != null) {
				databaseMolColumn = dName;
			} else {
				throw new InvalidSettingsException("No reference CDK compatible column in input table 2");
			}
		}

    	if ( refMolColumn == null
				|| (inSpecs[1].findColumnIndex(refMolColumn)) == -1) {
			String rName = null;
			for (DataColumnSpec s : inSpecs[1]) {
				if (s.getType().isCompatible(CDKValue.class) || s.getType().isCompatible(SmartsValue.class) ) {
					rName = s.getName();
				}
			}
			if (rName != null) {
				refMolColumn = rName;
			} else {
				throw new InvalidSettingsException("No reference CDK compatible column in input table 2");
			}
		}
    	
    	DataTableSpec outSpec = new DataTableSpec(createSpec(inSpecs[0]));
		return new DataTableSpec[] { new DataTableSpec(inSpecs[0], outSpec) };
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {

        settings.addInt(CFGKEY_TOPODISTANCECONSTRAINT, topologicalDistanceLimit);
        settings.addInt( CFGKEY_EXPTIMEOUT, expansionTimeOut );
        
    	settings.addString( CFGKEY_REFMOLCOLUMN, refMolColumn );  
    	settings.addString( CFGKEY_DATABASEMOLCOLUMN, databaseMolColumn );  
    	
    	settings.addString( CFGKEY_ALGORITHM, mappingAlgorithmName.toString() );  
    	settings.addString( CFGKEY_RETURNTYPE, returnType.toString() );  
    	settings.addString( CFGKEY_SIMILARITYTYPE, simType.toString() );  
    	settings.addString( CFGKEY_AGGREGATION, aggregationMethod.toString() );  
    	
    	settings.addBoolean( CFGKEY_OUTPUTMCS, outputMCS );
    	settings.addBoolean( CFGKEY_MOREGHOSTINFO, detailedGhostInfo );
    	settings.addBoolean( CFGKEY_RAYMONDHEURISTICS, useRaymondHeuristics );
    	settings.addBoolean( CFGKEY_RINGHEURISTICS, useRingHeuristics );
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // load (valid) settings from the config object.
        // It can be safely assumed that the settings are valided by the 
        // method below.
        
    	if( settings.containsKey( CFGKEY_REFMOLCOLUMN ) )
    		refMolColumn = settings.getString( CFGKEY_REFMOLCOLUMN );
    	
    	if( settings.containsKey( CFGKEY_DATABASEMOLCOLUMN ) )
    		databaseMolColumn = settings.getString( CFGKEY_DATABASEMOLCOLUMN );
    	
    	mappingAlgorithmName = ExtendedAlgorithm.valueOf( settings.getString( CFGKEY_ALGORITHM ) );
    	aggregationMethod = AggregationMethod.valueOf( settings.getString( CFGKEY_AGGREGATION ) );
    	returnType = ReturnType.valueOf( settings.getString( CFGKEY_RETURNTYPE ) );
    	simType = SimilarityType.valueOf( settings.getString( CFGKEY_SIMILARITYTYPE ) );
    	logger.warn( "aggregation = " + aggregationMethod );
    	System.out.println( "aggregation = " + aggregationMethod );
    	outputMCS = settings.getBoolean( CFGKEY_OUTPUTMCS );
    	detailedGhostInfo = settings.getBoolean( CFGKEY_MOREGHOSTINFO );
    	useRaymondHeuristics = settings.getBoolean(CFGKEY_RAYMONDHEURISTICS);
    	useRingHeuristics = settings.getBoolean(CFGKEY_RINGHEURISTICS);
    	topologicalDistanceLimit = settings.getInt(CFGKEY_TOPODISTANCECONSTRAINT);
    	expansionTimeOut = settings.getInt(CFGKEY_EXPTIMEOUT);
    	
    	if( outputMCS ) {
    		ghostBondIDsIndex++;
    		ghostSMARTSIndex++;
    	}
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
            
        // TODO check if the settings could be applied to our model
        // e.g. if the count is in a certain range (which is ensured by the
        // SettingsModel).
        // Do not actually set any values of any member variables.

        //m_count.validateSettings(settings);

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        
        // TODO load internal data. 
        // Everything handed to output ports is loaded automatically (data
        // returned by the execute method, models loaded in loadModelContent,
        // and user settings set through loadSettingsFrom - is all taken care 
        // of). Load here only the other internals that need to be restored
        // (e.g. data used by the views).

    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
       
        // TODO save internal models. 
        // Everything written to output ports is saved automatically (data
        // returned by the execute method, models saved in the saveModelContent,
        // and user settings saved through saveSettingsTo - is all taken care 
        // of). Save here only the other internals that need to be preserved
        // (e.g. data used by the views).

    }
    
    
    /**
	 * Provides a map of molecules and their corresponding rows.
	 * 
	 * @param bdt a buffered data table with DenseBitVector cells
	 * @param refColIndex a IAtomContainer (molecule) column index in the input table
	 * @return the map
	 */
	private Map<IAtomContainer, ArrayList<String>> getRefMols(DataTable dt, int refColIndex) {

		Map<IAtomContainer, ArrayList<String>> refs = new HashMap<IAtomContainer, ArrayList<String>>();
		SmartsParser sParser = new SmartsParser();
		
		for (DataRow row : dt) {
			// skip missing values
			if (row.getCell(refColIndex).isMissing()) {
				continue;
			}
			
			IAtomContainer refMol = null;
			if ( row.getCell( refColIndex ) instanceof SmartsValue ) {
				String refSMARTS =  ((SmartsValue) row.getCell( refColIndex )).getSmartsValue();
				//refMol =  SMARTSParser.parse( refSMARTS, DefaultChemObjectBuilder.getInstance() );
				
				refMol = sParser.parse( refSMARTS );
				ConvenienceTools.correctAtomBondTypes( (IQueryAtomContainer) refMol );
				//ConvenienceTools.countRings(refMol);
				ConvenienceTools.calculateImplicitHydrogens(refMol);
				//ConvenienceTools.initializeQuery( (IQueryAtomContainer) refMol );
				
				// XXX get rid of listeners - cause concurrency errors 
				refMol.setNotification(false);
				
				for( IBond b : refMol.bonds() )
					b.setNotification( false );
				
				for( IAtom a : refMol.atoms() )
					a.setNotification( false );
			
				
			//} else if( (((AdapterValue) row.getCell(refColIndex)).getAdapterError(CDKValue.class) != null) ) {  // CDKCell-compatible 
			} else if( row.getCell(refColIndex).getType().isCompatible( CDKValue.class ) ) {  // CDKCell-compatible 
				//CDKValue molValue = ((AdapterValue) row.getCell(refColIndex)).getAdapter(CDKValue.class);
				
				//refMol = ((CDKValue) row.getCell( refColIndex )).getAtomContainer();
				CDKValue molCell = ((AdapterValue) row.getCell(refColIndex)).getAdapter(CDKValue.class);
				refMol = molCell.getAtomContainer();
				
				logger.info("hs - " + refMol + " " + refMol.getAtomCount() + " " + refMol.getBondCount()  );
			} else {
				logger.error("ref molecule incompatible with CDK adapter values! ");
			}
			
			if (refs.containsKey(refMol)) {
				refs.get(refMol).add( row.getKey().getString() );  // store the row name
			} else {
				ArrayList<String> keyList = new ArrayList<String>();
				keyList.add(row.getKey().getString());
				refs.put(refMol, keyList);
			}
			
			
			// re-perceive aromaticity and ring information
			try {
				ConvenienceTools.countRings(refMol);  // essential for workings of modular product heuristics (Raymond-like ones)
				ConvenienceTools.calculateAromaticity( refMol );
			} catch (CDKException e) {
				e.printStackTrace();
			}
			
			refMol.setNotification(false);  // Hoping this will turn off listeners which cause concurrency problems
			
			// If bond origin information exists, output ghost substructure info in MCSs
			// currently error-prone - doesn't support mixtures of molecules (ones with and without the desired property)
			
	
			if( refMol.getBondCount() > 0 && refMol.getBond(0).getProperty(CDKSMARTSHyperstructureFitness.bondMolOriginType) != null ) {
				ghostSubstructureInfo = true;
				//System.out.println( "ghostSubstructureInfo is " + ghostSubstructureInfo );
			} else {
				detailedGhostInfo = false;
			}
			

		}
		return refs;
	}

	
	/**
	 * Provides a list of reference molecules in their given order.
	 * 
	 * @param bdt a buffered data table with DenseBitVector cells
	 * @param refColIndex a fingerprint column index in bdt
	 * @return the map
	 */
	@Deprecated
	private Map<IAtomContainer, ArrayList<String>> getMatrixRefs(DataTable dt, int refColIndex) {

		Map<IAtomContainer, ArrayList<String>> refs = new HashMap<IAtomContainer, ArrayList<String>>();
		SmartsParser sParser = new SmartsParser();
		
		for (DataRow row : dt) {
			/*if (row.getCell(refColIndex).isMissing()) {
				refs.add(null);
				continue;
			} */
			
			IAtomContainer refMol = null;
			if( row.getCell( refColIndex ) instanceof CDKValue ) {
				refMol = ((CDKValue) row.getCell( refColIndex )).getAtomContainer();
			} else if ( row.getCell( refColIndex ) instanceof SmartsValue ) {
				String refSMARTS =  ((SmartsValue) row.getCell( refColIndex )).getSmartsValue();
				//refMol = SMARTSParser.parse( refSMARTS, DefaultChemObjectBuilder.getInstance() );
				
				refMol = sParser.parse( refSMARTS );
				ConvenienceTools.correctAtomBondTypes( (IQueryAtomContainer) refMol );
				//ConvenienceTools.initializeQuery( (IQueryAtomContainer) refMol );
			}
			
			if (refs.containsKey(refMol)) {
				refs.get(refMol).add( row.getKey().getString() );  // store the row name
			} else {
				ArrayList<String> keyList = new ArrayList<String>();
				keyList.add(row.getKey().getString());
				refs.put(refMol, keyList);
			}
			
			//refs.add(refMol);
		}

		return refs;
	}

	/*
	private synchronized double calculateSimilarity( Isomorphism comparator, IAtomContainer rMol, IAtomContainer dbMol ) {
		try {
			comparator.init(rMol, dbMol, true, false);
		} catch (CDKException e) {
			
			e.printStackTrace();
		}
		
		double similarity = -1.0;
		try {
			similarity = comparator.getTanimotoBondSimilarity();
		} catch (IOException e) {
			
			e.printStackTrace();
		}
		
		return similarity;
	}
	
	private synchronized double calculateSimilarity( IAtomContainer rMol, IAtomContainer dbMol ) {
		
		Isomorphism comparator = null;
		if( mappingAlgorithmName == MappingAlgorithmName.VFLIB_cMCS )
			comparator = new Isomorphism(Algorithm.VFLibMCS, true);
		
		return calculateSimilarity( comparator, rMol, dbMol );
	}
	*/
    
    protected ExtendedCellFactory[] prepareExecute(final DataTable[] data) throws Exception {

    	// variable declaration
		final String sr = refMolColumn;
		startTime = System.currentTimeMillis();
		
		final int refColIndex = data[1].getDataTableSpec().findColumnIndex(sr);
		final Map<IAtomContainer, ArrayList<String>> mapRefs = getRefMols(data[1], refColIndex);
		rowCount = mapRefs.size();
		
		//rowCount = mapRefs.size();
		
		//logger.warn("There are " + rowCount + " reference molecules. "); 
		
		
		// choose algorithm to use
		/*
		if( mappingAlgorithmName == MappingAlgorithmName.VFLib_cMCS ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.VFLibMCS, true);
		} else if( mappingAlgorithmName == MappingAlgorithmName.CDKMCS ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.CDKMCS, true);
		} else if( mappingAlgorithmName == MappingAlgorithmName.vfLibGA_dMCS ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.vfLibGAdMCS, false);
		} else if( mappingAlgorithmName == MappingAlgorithmName.SubStructure ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.SubStructure, true);
		} else if( mappingAlgorithmName == MappingAlgorithmName.MCSPlus ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.MCSPlus, true);
		} else if( mappingAlgorithmName == MappingAlgorithmName.Zhu_AER_dMCS ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.ZhuAERdMCES, true);
		} else if( mappingAlgorithmName == MappingAlgorithmName.ChemAxon_cMCS ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.ChemAxon_cMCES, true);
		} else if( mappingAlgorithmName == MappingAlgorithmName.ChemAxon_dMCS ) {
			similarityComparator = new ExtendedIsomorphism(ExtendedAlgorithm.ChemAxon_dMCES, true);
		}
		*/
		
		final int databaseMolColIndex = data[0].getDataTableSpec().findColumnIndex( databaseMolColumn );
		final String weightProperty = "_weightSum";

		
		// assign weighted sums of bonds (if frequencies exist)
		for( IAtomContainer ref : mapRefs.keySet() ) {
				int weightSum = 0;
				
				if( ref == null )
					continue;
				
				// determine if there're weights or not (all bonds must have weights for this to work)
				useWeights = false;
				if( ref.getBondCount() > 0 && ref.getBond(0).getProperty(CDKSMARTSHyperstructureFitness.bondFrequencyType) != null )
					useWeights = true;
				
				/*useWeights = true;
				for( IBond b : ref.bonds() ) {
					if( b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ) == null ) {
						useWeights = false;
						break;
					}
				}*/
				
				
				if( useWeights ) {
					for( IBond b : ref.bonds() ) {
						//System.out.println( "freq: " + b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType ) );
						weightSum += (Integer) b.getProperty( CDKSMARTSHyperstructureFitness.bondFrequencyType );
					} 
				} else {
					weightSum = ref.getBondCount();
				}
				
				ref.setProperty( weightProperty, weightSum );
		}
		
		
		
		ExtendedCellFactory cf = new ExtendedCellFactory() {
			
			/*
			private double calculateSimilarity( ExtendedIsomorphism comparator, IAtomContainer dbMol, IAtomContainer rMol, boolean tversky ) {
				
				double similarity = -1.0;
				
				try {
					if( rMol instanceof IQueryAtomContainer ) {
						IQueryAtomContainer rMolQ = (IQueryAtomContainer) rMol;
						comparator.init(rMolQ, dbMol);
						logger.debug("IQueryAtomContainer used for similarity");
						System.out.println("IQueryAtomContainer used for similarity");
					} else {
						comparator.init(rMol, dbMol, true, false);
					}
					//comparator.setChemFilters(true, true, true);
					
					if( tversky )
						similarity = comparator.getTverskyBondSimilarity();
					else
						similarity = comparator.getTanimotoBondSimilarity();
					
				} catch (CDKException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e2) {
					// TODO Auto-generated catch block
					e2.printStackTrace();
				} catch (NullPointerException e3) {
					e3.printStackTrace();
					logger.warn("Null pointer - similarity set to 0 " + rMol.getBondCount() + " " + dbMol.getBondCount() + " " + comparator.getAllBondMaps());
					similarity = 0.0;
				}
				
				return similarity;
			}*/
			
			private void colourMCSAtoms( Map<IAtom, IAtom> mapping, IAtomContainer dbMol ) {
				
				Color[] color = CDKNodeUtils.generateColors(1);
				
				for( IAtom at : mapping.keySet() ) {
					at.setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
				}
				
				for( IAtom at : mapping.values() ) {  // dbmol atoms
					at.setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
					//System.out.println( key +  " bondt: " + dbMol.contains(key) + " " + dbMol.contains(key.getAtom(0)) + " " + dbMol.contains(key.getAtom(1)) + " " + rMol.contains(key.getAtom(0)) + " " + key.getProperties() );
					
				}
				
				dbMol.getAtom(0).setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
				dbMol.getAtom(1).setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
				dbMol.getAtom(4).setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
			}
			
			
			private void colourMCSAtoms( Collection<Integer> indices, IAtomContainer mol ) {
				
				Color[] color = CDKNodeUtils.generateColors(1);
				//System.out.println(" starting atom map ");
				 
				for( Integer index : indices ) {  // dbmol atoms
					
					IAtom at = mol.getAtom( index );
					at.setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
					
					at.setProperty(StandardGenerator.HIGHLIGHT_COLOR, color[0]  );
					
					//System.out.println(  at +  " atom: " + mol.contains(at) + " " + at.getProperties() );
					
				}
				
				//dbMol.getAtom(0).setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());  // a test
				
			}
			
			
			private void colourMCSBonds( Map<IBond, IBond> mapping, IAtomContainer dbMol ) {
				
				Color[] color = CDKNodeUtils.generateColors(1);
				Color red = Color.RED;
				
		
				
				
				for( IBond bond : mapping.keySet() ) {
					bond.setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
					//key.getAtom(0).setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
					//System.out.println( bond +  " bondt: " + dbMol.contains(bond) + " " + dbMol.contains(bond.getAtom(0)) + " " + dbMol.contains(bond.getAtom(1)) + " " + bond.getProperties() );
					
				}
			 
				
				/*dbMol.getBond(0).setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
				dbMol.getBond(2).setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
				dbMol.getBond(3).setProperty(CDKConstants.ANNOTATIONS, red.getRGB());*/
				
				//System.out.println("bond map finished");
				/*
				for( IBond value : mapping.values() ) {
					IBond b = rMol.getBond(value.getAtom(0), value.getAtom(1));
					b.setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
				}
				*/
			}
			
			
			private void colourMCSBonds( Collection<Integer> indices, IAtomContainer mol ) {
				
				Color[] color = CDKNodeUtils.generateColors(1);
				//System.out.println(" starting bond map ");
				 
				for( Integer index : indices ) {  // dbmol atoms
					
					IBond b = mol.getBond( index );
					b.setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
					//e.getValue().setProperty(StandardGenerator.HIGHLIGHT_COLOR, color[i]);
					b.setProperty(StandardGenerator.HIGHLIGHT_COLOR, color[0]  );
					
					//System.out.println(  b +  " bond: " + mol.contains(b) + " " + b.getProperties() );
					
				}
				
			}
			
			/*
			private void colourMCSBonds( Map<Integer, Integer> mapping, IAtomContainer dbMol ) {
				
				Color[] color = CDKNodeUtils.generateColors(2);
				Color red = Color.RED;
				// generate table to translate String IDs to Bonds, seeing as the Mapping process destroys the relation with the original molecule
				
				
				for( Entry<Integer, Integer> pair : mapping.entrySet() ) {
					
					//IBond hsBond = rMol.getBond( pair.getKey() );
					IBond dbBond = dbMol.getBond( pair.getValue() );
					
					//hsBond.setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
					dbBond.setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
					//System.out.println( dbBond +  " bondt: " + dbMol.contains(dbBond.getAtom(0)) + " " + dbMol.contains(dbBond.getAtom(1)) + " " + rMol.contains(hsBond.getAtom(0)) );
				}
				
				dbMol.getAtom(0).setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
				dbMol.getAtom(6).setProperty(CDKConstants.ANNOTATIONS, red.getRGB());
				System.out.println("bond map finished");
				
				for( IBond value : mapping.values() ) {
					IBond b = rMol.getBond(value.getAtom(0), value.getAtom(1));
					b.setProperty(CDKConstants.ANNOTATIONS, color[0].getRGB());
				}
				
			}
			*/
			
			
			@Override
			public DataCell[] getCells(final DataRow row) {
				
				ExtendedIsomorphism sComparator = null;
				MCSMethods mcsMapper = null;
				
				
				
			
				
				

				DataCell molCell = row.getCell(databaseMolColIndex);
				DataCell[] cells = new DataCell[getColumnSpecs().length];
				if (molCell.isMissing()) {
					Arrays.fill(cells, DataType.getMissingCell());
					return cells;
				}
				if (!(molCell instanceof CDKValue)) {
					throw new IllegalArgumentException("No Database Molecule cell at " + databaseMolColIndex + ": "
							+ molCell.getClass().getName());
				}
				
				IAtomContainer dbMol = null;
				
				
				
				// Cloning is performed so that the original molecules aren't modified
				try {
					//dbMol = CDKNodeUtils.getExplicitClone( ((CDKValue) row.getCell( databaseMolColIndex )).getAtomContainer() );  // don't use as it adds hydrogens
					dbMol = ((CDKValue) molCell).getAtomContainer().clone();
				} catch (CloneNotSupportedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				// re-perceive aromaticity and ring information
				try {
					AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(dbMol);
					ConvenienceTools.calculateImplicitHydrogens(dbMol);
					ConvenienceTools.correctAtomBondTypes(dbMol);
					ConvenienceTools.countRings( dbMol );
					ConvenienceTools.calculateAromaticity( dbMol );
				} catch (CDKException e) {
					e.printStackTrace();
				}
				
				
				/*sComparator = new ExtendedIsomorphism(mappingAlgorithmName, true);
				//sComparator = new ExtendedIsomorphism(refMol, dbMol, mappingAlgorithmName, true, false, true, useRaymondHeuristics, topologicalDistanceLimit, expansionTimeOut);
				sComparator.setTopologicalDistanceLimit(topologicalDistanceLimit);
				sComparator.setUseRaymondHeuristics(useRaymondHeuristics);
				sComparator.setTimeLimit( expansionTimeOut );*/
				//sComparator = new ExtendedIsomorphism(ExtendedAlgorithm.BK, true);
				
				SimilarityComparator simHub = new SimilarityComparator( sComparator, mcsMapper, useWeights, ghostSubstructureInfo );
				
				
				
				ArrayList<String> pkey = null;
				
				if ( aggregationMethod == AggregationMethod.Matrix) {
					List<DataCell> resultsSimilarity = new ArrayList<DataCell>();
					List<DataCell> resultsTanimoto = new ArrayList<DataCell>();
					List<DataCell> resultsNumBonds = new ArrayList<DataCell>();
					List<DataCell> resultsNumRefBonds = new ArrayList<DataCell>();
					List<DataCell> resultsNumRefWeights = new ArrayList<DataCell>();
					IntCell resultsNumDbBonds = new IntCell(0);
					List<DataCell> resultsMCSWeights = new ArrayList<DataCell>();
					List<DataCell> resultsMCSTopologies = new ArrayList<DataCell>();
					List<DataCell> resultsFragmentSizes = new ArrayList<DataCell>();
					List<DataCell> resultsMCSSMARTS = new ArrayList<DataCell>();
					List<DataCell> resultsMCSTimes = new ArrayList<DataCell>();
					List<DataCell> resultsModProdTimes = new ArrayList<DataCell>();
					List<DataCell> resultsModProdSizes = new ArrayList<DataCell>();
					List<DataCell> resultsModProdDensities = new ArrayList<DataCell>();
					List<DataCell> resultsMCSUniqueSubstructures = new ArrayList<DataCell>();
					List<DataCell> resultsMCSUniqueGhostSubstructures = new ArrayList<DataCell>();
					List<DataCell> resultsMCSGhostSMARTSList = new ArrayList<DataCell>();
					
					
					pkey = new ArrayList<String>();
					Iterator<Map.Entry<IAtomContainer, ArrayList<String>>> it = mapRefs.entrySet().iterator();
					
					while (it.hasNext()) {
						
						Map.Entry<IAtomContainer, ArrayList<String>> pairs = it.next();
						IAtomContainer refMol = pairs.getKey();
						
						if (refMol == null ||  dbMol == null || dbMol.getAtomCount() < 1 ) { 
							resultsSimilarity.add(DataType.getMissingCell());
							resultsTanimoto.add(DataType.getMissingCell());
							resultsNumBonds.add(DataType.getMissingCell());
							resultsNumRefBonds.add(DataType.getMissingCell());
							resultsMCSWeights.add( DataType.getMissingCell() );
							resultsMCSTopologies.add( DataType.getMissingCell() );
							resultsFragmentSizes.add(DataType.getMissingCell()); 
							resultsMCSTimes.add(DataType.getMissingCell()); 
							resultsModProdTimes.add(DataType.getMissingCell()); 
							resultsModProdSizes.add(DataType.getMissingCell()); 
							resultsModProdDensities.add(DataType.getMissingCell()); 
							resultsMCSUniqueSubstructures.add( DataType.getMissingCell() );
							resultsMCSUniqueGhostSubstructures.add( DataType.getMissingCell() );
							resultsMCSGhostSMARTSList.add( DataType.getMissingCell() );
							
							logger.warn( "calculateSimilarity attempted here and failed: " + row.getCell(0) + " against reference " + pairs.getValue()  );
						} else {
							//sComparator.setChemFilters(false, false, false);
							simHub.calculateSimilarity(refMol, dbMol );
							
							
							double similarity = simHub.getSimilarityValue( simType ) ;
							double tanimoto = simHub.tanimoto;
							
							
							if( simHub.bondMaps != null && simHub.bondMaps.size() > 0  ) {
							
								//ArrayList<Integer> fSizes = new ArrayList<Integer>( sComparator.getFragmentSizes().length );
								ArrayList<DataCell> fSizes = new ArrayList<DataCell>( simHub.fragmentSizes.length );
								for( int i : simHub.fragmentSizes )  fSizes.add( new IntCell( i ) );
								
								ArrayList<DataCell> mcsWeights = new ArrayList<DataCell>( simHub.commonBondWeights.length );
								for( int i : simHub.commonBondWeights )  mcsWeights.add( new IntCell( i ) );
								
								ArrayList<DataCell> mcsTopologies = new ArrayList<DataCell>( simHub.commonBondTopologies.length );
								for( String i : simHub.commonBondTopologies )  mcsTopologies.add( new StringCell( i ) );
								
								
								
								int commonBonds = simHub.bondMaps.get(0).size();
								
								resultsSimilarity.add(new DoubleCell( similarity ));
								resultsTanimoto.add(new DoubleCell( tanimoto ));
								resultsNumBonds.add(new IntCell( commonBonds ));
								resultsNumRefBonds.add(new IntCell( refMol.getBondCount() ) ) ;
								resultsNumRefWeights.add( new IntCell( (Integer) refMol.getProperty(weightProperty) ) );
								
								resultsNumDbBonds = new IntCell( dbMol.getBondCount() ) ;
								resultsMCSWeights.add( CollectionCellFactory.createListCell( mcsWeights ) );
								resultsMCSTopologies.add( CollectionCellFactory.createListCell( mcsTopologies ) );
								resultsFragmentSizes.add( CollectionCellFactory.createListCell( fSizes ) );
								
								resultsMCSTimes.add( new IntCell( (int) simHub.mcsExecTime ) );
								resultsModProdTimes.add( new IntCell( (int) simHub.modProdConstructionTime ) ); 
								resultsModProdSizes.add( new IntCell( (int) simHub.modProdNodeCount ) ); 
								resultsModProdDensities.add( new DoubleCell( simHub.modProdEdgeDensity ) ); 
								
								resultsMCSUniqueSubstructures.add( new IntCell( simHub.uniqueSubstructureCount ) );
								resultsMCSUniqueGhostSubstructures.add( new IntCell( simHub.mcsGhostCount ) );
								
								
								if( detailedGhostInfo && simHub.mcsGhostSMARTSSet.size() > 0 ) {
									StringBuilder ghostSMARTSSet = new StringBuilder( simHub.mcsGhostSMARTSSet.size() * simHub.mcsGhostSMARTSSet.get(0).length() );
									ghostSMARTSSet.append( simHub.mcsGhostSMARTSSet.get(0) );
									
									for( int i = 1; i < simHub.mcsGhostSMARTSSet.size(); i++ ) {
										ghostSMARTSSet.append(".");
										ghostSMARTSSet.append( simHub.mcsGhostSMARTSSet.get(i) );
									}
									
									resultsMCSGhostSMARTSList.add( new SmartsCell( ghostSMARTSSet.toString() ) );
								}
								
								
								
								if( outputMCS ) {
										resultsMCSSMARTS.add( new SmartsCell( simHub.mcsSMARTS ) );
								}
								//resultsFragmentSizes.add(DataType.getMissingCell());
								
								//logger.warn( "calculateSimilarity called here " + sComparator.isSubgraph() + " " + commonBonds + " " + dbMol.getAtomCount() + " " + cells.length  );
							} else {
								resultsSimilarity.add(DataType.getMissingCell());
								resultsTanimoto.add(DataType.getMissingCell());
								resultsNumBonds.add(DataType.getMissingCell());
								resultsNumRefBonds.add(DataType.getMissingCell());
								resultsFragmentSizes.add(DataType.getMissingCell());
							}
							pkey.add( pairs.getValue().get(0) );
						}
					}
					
					cells[tverskyCIndex] = CollectionCellFactory.createListCell(resultsSimilarity);
					cells[tanimotoCIndex] = CollectionCellFactory.createListCell(resultsTanimoto);
					cells[commonBondsCIndex] = CollectionCellFactory.createListCell(resultsNumBonds);
					cells[refBondsCIndex] = CollectionCellFactory.createListCell(resultsNumRefBonds);
					cells[refWeightSumCIndex] = CollectionCellFactory.createListCell(resultsNumRefWeights);
					cells[dbBondsCIndex] = resultsNumDbBonds;
					cells[mcsWeightsCIndex] = CollectionCellFactory.createListCell(resultsMCSWeights);
					cells[mcsTopologiesCIndex] = CollectionCellFactory.createListCell(resultsMCSTopologies);
					cells[fragmentSizesCIndex] = CollectionCellFactory.createListCell(resultsFragmentSizes);
					cells[mcsTimeCIndex] = CollectionCellFactory.createListCell(resultsMCSTimes);
					cells[modProdTimeCIndex] = CollectionCellFactory.createListCell(resultsModProdTimes);
					cells[modProdSizeCIndex] = CollectionCellFactory.createListCell(resultsModProdSizes);
					cells[modProdDensityCIndex] = CollectionCellFactory.createListCell(resultsModProdDensities);
					cells[uniqueSubstructuresCIndex] = CollectionCellFactory.createListCell(resultsMCSUniqueSubstructures);
					cells[uniqueGhostSubstructuresCIndex] = CollectionCellFactory.createListCell(resultsMCSUniqueGhostSubstructures);
				
					
					if( outputMCS ) {
						//if( mappingAlgorithmName == MappingAlgorithmName.ChemAxon_dMCS || mappingAlgorithmName == MappingAlgorithmName.ChemAxon_dMCS ) {
							cells[outputMCSIndex] = CollectionCellFactory.createListCell( resultsMCSSMARTS );
						//} else {
						//	cells[outputMCSIndex] = DataType.getMissingCell();
						//}
					}
					
					if( detailedGhostInfo ) { 
						cells[ ghostSMARTSIndex ] = CollectionCellFactory.createListCell( resultsMCSGhostSMARTSList );
					}
					
				} else {
					double coeff = 0.0f;
					double pcoeff = -0.05f;
					double tanCoeff = -0.05f;
					int commonBonds = 0, refBonds = 0, dbBonds = 0;
					int uniqueSubstructures = 0, uniqueGhostSubstructures = 0;
					
					Iterator<Map.Entry<IAtomContainer, ArrayList<String>>> it = mapRefs.entrySet().iterator();
					ArrayList<DataCell> fSizes = null;
					List<DataCell> mcsWeights = null; //= new ArrayList<DataCell>();
					List<DataCell> mcsTopologies = null;
					StringBuilder resultsMCSGhostSMARTSList =  new StringBuilder( 100 );
					int refWeightSum = -1;
					
					int mcsTime = 0;
					int modProdTime = 0;
					int modProdSize = 0;
					double modProdDensity = 0.0;
					 
					Map<IBond, IBond> bestMap = null;
					//AtomAtomMapping bestAtomMap = null;
					Map<Integer, Integer> bestAtomIndexMap = null;
					List<int[]> bestBondIndexMap = null;
					IAtomContainer bestRef = null;
					
					if ( aggregationMethod == AggregationMethod.Minimum) {
						
						pcoeff = 1;
						
						
						while (it.hasNext()) {
							
							Map.Entry<IAtomContainer, ArrayList<String>> pairs = it.next();
							IAtomContainer refMol = pairs.getKey();
							/*
							try {
								db2 = CDKNodeUtils.getExplicitClone(dbMol);
								refMol = CDKNodeUtils.getExplicitClone(refMol);
							} catch (CDKException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
							*/
							//coeff = Tanimoto.calculate(bs, pairs.getKey());
							
							
							// ExtendedIsomorphism sC2 = sComparator;
							simHub.calculateSimilarity( refMol, dbMol );
							//coeff = calculateSimilarity( sC2, dbMol, refMol, true );
							coeff = simHub.getSimilarityValue( simType ) ;
							
							
							if( dbMol == null || dbMol.getAtomCount() < 1 || simHub.bondMaps == null || simHub.bondMaps.size() == 0 ) {
								logger.warn( "Min calculateSimilarity attempted here and failed: " + row.getCell(0) + " against reference " + pairs.getValue() );
								continue;
							}
							
							//if( sC2.getFirstMapping() != null )
							//	logger.warn( "Min calculateSimilarity called here " + sC2.isSubgraph() + " " + sC2.getFirstMapping().size() + " " + dbMol.getAtomCount() + " " + cells.length + " " + row.getCell(0) /*+ sComparator.getAllBondMaps().get(0).size()*/  );
							
							
							dbBonds = dbMol.getBondCount();
							
							
							if (coeff <= pcoeff && simHub.bondMaps != null ) {
								refBonds = refMol.getBondCount();
								bestRef = refMol;
								
								pcoeff = coeff;
								pkey = pairs.getValue();
								bestMap = new HashMap<IBond, IBond>( simHub.bondMaps.get(0) );  // clone the map as CDK often destroys the original
								//bestAtomMap = simHub.atomMaps.get(0);
								bestAtomIndexMap = simHub.atomIndexMaps.get(0);
								bestBondIndexMap = simHub.bondIndexMaps.get(0);
								//System.out.println( "atom indices - " + bestAtomIndexMap );
								
								tanCoeff = simHub.tanimoto;
								commonBonds = bestMap.size();
								
								refWeightSum = bestRef.getProperty( weightProperty );
								
								
								fSizes = new ArrayList<DataCell>( simHub.fragmentSizes.length );
								for( int i : simHub.fragmentSizes )  fSizes.add( new IntCell( i ) );
								
								mcsWeights = new ArrayList<DataCell>( simHub.commonBondWeights.length );
								for( int i : simHub.commonBondWeights )  mcsWeights.add( new IntCell( i ) );
								
								mcsTopologies = new ArrayList<DataCell>( simHub.commonBondTopologies.length );
								for( String i : simHub.commonBondTopologies )  mcsTopologies.add( new StringCell( i ) );
								
								mcsTime = (int) simHub.mcsExecTime;
								modProdTime = (int) simHub.modProdConstructionTime;
								modProdSize = simHub.modProdNodeCount;
								modProdDensity = simHub.modProdEdgeDensity;
								
								uniqueSubstructures = simHub.uniqueSubstructureCount;
								uniqueGhostSubstructures = simHub.mcsGhostCount;
								
								if( detailedGhostInfo && simHub.mcsGhostSMARTSSet.size() > 0 ) {
									resultsMCSGhostSMARTSList.append( simHub.mcsGhostSMARTSSet.get(0) );
									
									for( int i = 1; i < simHub.mcsGhostSMARTSSet.size(); i++ ) {
										resultsMCSGhostSMARTSList.append(".");
										resultsMCSGhostSMARTSList.append( simHub.mcsGhostSMARTSSet.get(i) );
									}
								}
								
								
								//logger.warn( "calculateSimilarity called here " + similarityComparator.isSubgraph() + " " + similarityComparator.getFirstMapping().size() + " " + dbMol.getAtomCount() + " " + cells.length  );
							}
						}
						
						//logger.warn( "calculateSimilarity called here " + similarityComparator.isSubgraph() + " " + bestMap.size() + " " + dbMol.getAtomCount()  );
						
					} else if ( aggregationMethod == AggregationMethod.Maximum) {
						
						while (it.hasNext()) {
							Map.Entry<IAtomContainer, ArrayList<String>> pairs = it.next();
							IAtomContainer refMol = pairs.getKey();
							
							simHub.calculateSimilarity( refMol, dbMol );
							coeff = simHub.getSimilarityValue( simType ) ;
							
							if( dbMol == null || dbMol.getAtomCount() < 1 || simHub.bondMaps == null || simHub.bondMaps.size() == 0 ) {
								logger.warn( "Min calculateSimilarity attempted here and failed: " + row.getCell(0) + " against reference " + pairs.getValue()   );
								continue;
							}
							
							
							
							dbBonds = dbMol.getBondCount();
							
							if (coeff >= pcoeff && simHub.bondMaps != null ) {
								refBonds = refMol.getBondCount();
								bestRef = refMol;
								
								pcoeff = coeff;
								pkey = pairs.getValue();
								bestMap = new HashMap<IBond, IBond>( simHub.bondMaps.get(0) );  // clone the map as CDK often destroys the original
								//bestAtomMap = simHub.atomMaps.get(0);
								bestAtomIndexMap = simHub.atomIndexMaps.get(0);
								bestBondIndexMap = simHub.bondIndexMaps.get(0);
								
								tanCoeff = simHub.tanimoto;
								commonBonds = bestMap.size();
								
								refWeightSum = bestRef.getProperty( weightProperty );
								
								
								fSizes = new ArrayList<DataCell>( simHub.fragmentSizes.length );
								for( int i : simHub.fragmentSizes )  fSizes.add( new IntCell( i ) );
								
								mcsWeights = new ArrayList<DataCell>( simHub.commonBondWeights.length );
								for( int i : simHub.commonBondWeights )  mcsWeights.add( new IntCell( i ) );
								
								mcsTopologies = new ArrayList<DataCell>( simHub.commonBondTopologies.length );
								for( String i : simHub.commonBondTopologies )  mcsTopologies.add( new StringCell( i ) );
								
								mcsTime = (int) simHub.mcsExecTime;
								modProdTime = (int) simHub.modProdConstructionTime;
								modProdSize = simHub.modProdNodeCount;
								modProdDensity = simHub.modProdEdgeDensity;
								
								uniqueSubstructures = simHub.uniqueSubstructureCount;
								uniqueGhostSubstructures = simHub.mcsGhostCount;
								
								if( detailedGhostInfo && simHub.mcsGhostSMARTSSet.size() > 0 ) {
									resultsMCSGhostSMARTSList.append( simHub.mcsGhostSMARTSSet.get(0) );
									
									for( int i = 1; i < simHub.mcsGhostSMARTSSet.size(); i++ ) {
										resultsMCSGhostSMARTSList.append(".");
										resultsMCSGhostSMARTSList.append( simHub.mcsGhostSMARTSSet.get(i) );
									}
								}
							}
						}
						 
						//try {
							//IQueryAtomContainer dbm = QueryAtomContainerCreator.createAnyAtomAnyBondContainer(dbMol, false);
							//sC2.init(dbMol, bestRef, true, false);
							
							//logger.warn( "calculateSimilarity called here " + sComparator.isSubgraph() + " " + sComparator.getAllAtomIndexMapping().get(0) + " " + dbMol.getAtomCount() + " " + cells.length + " " + bestMap  );
							//logger.warn(  sgTester.isSubgraph( dbMol, QueryAtomContainerCreator.createAnyAtomAnyBondContainer( bestRef, true ) ) );
						//} catch (CDKException e) {
							// TODO Auto-generated catch block
						//	e.printStackTrace();
						//}
						
							
						
					} else if ( aggregationMethod == AggregationMethod.Average ) {
						tanCoeff = 0;
						commonBonds = 0;
						refBonds = 0;
						dbBonds = 0;
						
						fSizes = new ArrayList<DataCell>(0);
						
						mcsWeights = new ArrayList<DataCell>(  );
						//for( int i : simHub.commonBondWeights )  mcsWeights.add( new IntCell( i ) );
						
						mcsTopologies = new ArrayList<DataCell>( );
						
						while (it.hasNext()) {
							
							Map.Entry<IAtomContainer, ArrayList<String>> pairs = it.next();
							//coeff += calculateSimilarity( sComparator, dbMol, pairs.getKey(), true );
							
							simHub.calculateSimilarity( pairs.getKey(), dbMol );
							tanCoeff += simHub.getSimilarityValue( simType ) ;
								
							//tanCoeff += ConvenienceTools.calculateTversky( sComparator.getAllBondMaps().get(0).size(), pairs.getKey().getBondCount(), dbMol.getBondCount(), 1.0, 1.0 );
							
							
							//commonBonds += (mappingAlgorithmName == MappingAlgorithmName.ChemAxon_dMCS || mappingAlgorithmName == MappingAlgorithmName.ChemAxon_cMCS) ? sComparator.getCAMCSSize() : sComparator.getFirstBondMap().size();
							bestMap = new HashMap<IBond, IBond>( simHub.bondMaps.get(0) );  // clone the map as CDK often destroys the original
							commonBonds += bestMap.size();
							
							//commonBonds += similarityComparator.getFirstBondMap().size();
							refBonds += pairs.getKey().getBondCount();
							dbBonds += dbMol.getBondCount();
							
							mcsTime += (int) simHub.mcsExecTime;
							modProdTime += (int) simHub.modProdConstructionTime;
							modProdSize += simHub.modProdNodeCount;
							modProdDensity += simHub.modProdEdgeDensity;
							
							uniqueSubstructures += simHub.uniqueSubstructureCount;
							uniqueGhostSubstructures += simHub.mcsGhostCount;
								
								//ArrayList<DataCell> fSizesElement = new ArrayList<DataCell>( sComparator.getFragmentSizes().length );
								//for( int i : sComparator.getFragmentSizes() )  fSizesElement.add( new IntCell( i ) );
								
								//fSizes.add( CollectionCellFactory.createListCell( fSizesElement ) );
							
						}
						pcoeff = coeff / rowCount;
						tanCoeff /= rowCount;
						commonBonds /= rowCount;
						refBonds /= rowCount;
						dbBonds /= rowCount;
						
						mcsTime /= rowCount;
						modProdTime /= rowCount;
						modProdSize /= rowCount;
						modProdDensity /= rowCount;
						
						uniqueSubstructures /= rowCount;
						uniqueGhostSubstructures /= rowCount;
						
						pkey = new ArrayList<String>();
						
					}
					
					
					//logger.warn("mcsWeights " + mcsWeights + " " + simHub.commonBondWeights.length  );
					cells[tverskyCIndex] = new DoubleCell(pcoeff);
					cells[tanimotoCIndex] = new DoubleCell(tanCoeff);
					cells[commonBondsCIndex] = new IntCell( commonBonds );
					cells[refBondsCIndex] = new IntCell( refBonds );
					cells[refWeightSumCIndex] = new IntCell( refWeightSum );
					cells[dbBondsCIndex] = new IntCell( dbBonds );
					
					if( commonBonds > 0 ) {
						try {
							cells[mcsWeightsCIndex] = CollectionCellFactory.createListCell( mcsWeights );
							cells[mcsTopologiesCIndex] = CollectionCellFactory.createListCell( mcsTopologies );
							
							cells[mcsTimeCIndex] = new IntCell( mcsTime );
							cells[modProdTimeCIndex] = new IntCell( modProdTime );
							cells[modProdSizeCIndex] = new IntCell( modProdSize );
							cells[modProdDensityCIndex] = new DoubleCell( modProdDensity );
							
							cells[uniqueSubstructuresCIndex] = new IntCell( uniqueSubstructures );
							cells[uniqueGhostSubstructuresCIndex] = new IntCell( uniqueGhostSubstructures );
						} catch (Exception e1) {
							// TODO Auto-generated catch block
							logger.error("Error - mcsWeights is null pointer at row " + row.getCell(0) );
							e1.printStackTrace();
						}
						
						try {
							cells[fragmentSizesCIndex] = CollectionCellFactory.createListCell( fSizes );
						} catch (NullPointerException e) {
							System.out.println("for some reason fSizes is " + fSizes );
							logger.warn( "for some reason fSizes is " + fSizes + " | " + row.getCell(0) + simHub.comparator.isSubgraph() + " " + simHub.comparator.getAllBondIndexMaps() + " " + dbMol.getAtomCount() + " " + cells.length );
							cells[fragmentSizesCIndex] = CollectionCellFactory.createListCell( new ArrayList<DataCell>(0) );
						} 
						
						if( detailedGhostInfo ) { 
							cells[ ghostSMARTSIndex ] = new SmartsCell( resultsMCSGhostSMARTSList.toString() );
						}
					} else {  // empty cells
						cells[mcsWeightsCIndex] = CollectionCellFactory.createListCell( new ArrayList<DataCell>(0) );
						cells[mcsTopologiesCIndex] = CollectionCellFactory.createListCell( new ArrayList<DataCell>(0) );
						cells[fragmentSizesCIndex] = CollectionCellFactory.createListCell( new ArrayList<DataCell>(0) );
						cells[mcsTimeCIndex] = DataType.getMissingCell();
						cells[modProdTimeCIndex] = DataType.getMissingCell();
						cells[modProdSizeCIndex] = DataType.getMissingCell();
						cells[modProdDensityCIndex] = DataType.getMissingCell();
						cells[uniqueSubstructuresCIndex] = DataType.getMissingCell();
						cells[uniqueGhostSubstructuresCIndex] = DataType.getMissingCell();
						
						if( detailedGhostInfo ) { 
							cells[ ghostSMARTSIndex ] = DataType.getMissingCell();
						}
					}
					
					
					
					if( outputMCS ) {
						if( bestBondIndexMap != null  && bestBondIndexMap.size() > 0 ) {
								//colourMCSBonds(bestMap, dbMol);  // Does this actually colour bonds?
								//colourMCSAtoms(bestAtomMap, dbMol);
								
								IAtomContainer refCopy = null;
								try {
									refCopy = bestRef.clone();
								} catch (CloneNotSupportedException e) {
									e.printStackTrace();
								}
								
								for( IAtom at : refCopy.atoms() ) {
									at.setProperty(CDKConstants.ANNOTATIONS, null);
								}
								
								for( IBond b : refCopy.bonds() ) {
									b.setProperty(CDKConstants.ANNOTATIONS, null);
								}
								
								// XXX  As of CDK 1.5.3 seems to colour both atoms and bonds, despite no bonds being specified
								// Hence the use of the CDKCell3 class for rendering
								colourMCSAtoms( bestAtomIndexMap.keySet(), refCopy ); 
								colourMCSAtoms( bestAtomIndexMap.values(), dbMol   ); 
								
								List<Integer> refMolBondIndices = new ArrayList<Integer>( bestBondIndexMap.size() );
								List<Integer> dbMolBondIndices = new ArrayList<Integer>( bestBondIndexMap.size() );
								for( int[] pair : bestBondIndexMap ) {
									refMolBondIndices.add( pair[0] );
									dbMolBondIndices.add( pair[1] );
								}
								
								colourMCSBonds( refMolBondIndices, refCopy); 
								colourMCSBonds( dbMolBondIndices, dbMol); 
								
								//dbMol.getBond(2).setProperty(CDKConstants.ANNOTATIONS, 65535);
								dbMol.add(refCopy);
								cells[outputMCSIndex] = CDKCell3.createCDKCell(dbMol);
						} else {
							cells[outputMCSIndex] = DataType.getMissingCell();
						}
					}
				}
				
				
				
				
				// add reference molecule names
				List<StringCell> res = new ArrayList<StringCell>();
				if( pkey != null && pkey.size() > 0 ) {
					for (String st : pkey ) {
						res.add(new StringCell(st));
					}
				}

				// as above
				if (res.size() > 0) {
					if ( returnType.equals(ReturnType.String)) {
						if (res.size() == 1)
							cells[refNamesCIndex] = res.get(0);
						else {
							String resString = "";
							for (StringCell cell : res) {
								resString += (cell.getStringValue() + "|");
							}
							resString = resString.substring(0, resString.lastIndexOf("|"));
							cells[refNamesCIndex] = new StringCell(resString);
						}
					} else if ( returnType.equals(ReturnType.Collection)) {
						cells[refNamesCIndex] = CollectionCellFactory.createListCell(res);
					}
				} else {
					cells[refNamesCIndex] = DataType.getMissingCell();
				}
				
				
				/*IAtomContainer dbMolOrig = ((CDKValue) row.getCell( databaseMolColIndex )).getAtomContainer();
				dbMolOrig.getBond(2).setProperty(CDKConstants.ANNOTATIONS, -65536);
				*/
				/*for( IBond bond : dbMol.bonds() ) {
					System.out.println( bond.getProperties() );
				}
				
				for( IAtom atom : dbMol.atoms() ) {
					System.out.println( "atom - " + atom.getProperties() );
				}*/
				
				//logger.info(cells[0] + " " + cells[11]  );
				return cells;
			}

			
			
			@Override
			public ColumnDestination[] getColumnDestinations() {

				return new ColumnDestination[] { new AppendColumn() };
			}

			@Override
			public DataColumnSpec[] getColumnSpecs() {

				return createSpec(data[0].getDataTableSpec());
			}
		};

		
		
		return new ExtendedCellFactory[] { cf };
	}


	@Override
	protected BufferedDataTable[] postExecute(BufferedDataTable[] res,
			ExecutionContext exec) {
		
		BufferedDataTable[] temp = super.postExecute(res, exec);
		
		long execTime = System.currentTimeMillis() - startTime;
		pushFlowVariableInt("execTime", (int) execTime);
		
		return temp;
	}
    
	
	
    private class SimilarityComparator {
    	
    	public SimilarityComparator( ExtendedIsomorphism e, MCSMethods m, boolean weights, boolean ghostSubstructures ) {
    		this.comparator = e;
    		this.mcsMapper = m;
    		this.useWeights = weights;
    		this.ghostSubstructures = ghostSubstructures;
    		
    		this.uniqueSubstructureCount = 0;
    		this.mcsGhostCount = 0;
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
					
					
					comparator = new ExtendedIsomorphism(
							rMol, 
							dbMol, 
							mappingAlgorithmName, 
							true, 
							false, 
							true, 
							useRaymondHeuristics,
							useRingHeuristics,
							topologicalDistanceLimit, 
							expansionTimeOut);
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
					logger.warn("Null pointer - similarity set to -1.0 " + rMol.getBondCount() + " " + dbMol.getBondCount() );
					logger.warn( "common bonds - " + comparator.getAllBondMaps() );
					logger.warn("settings - " + mappingAlgorithmName  + " " + useRaymondHeuristics + " " + topologicalDistanceLimit + " " + expansionTimeOut  );
					logger.warn("stats - " + fragmentSizes  + " " + mcsSMARTS + " " + modProdConstructionTime + " " + modProdNodeCount + " " + modProdEdgeDensity + " " + mcsExecTime  );
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
	    			
	    			
	    			if( detailedGhostInfo )
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
	    					
	    					if( detailedGhostInfo && isGhost && ! uniqueSubstructures.containsKey(identifier) )
	    						if( test instanceof IQueryAtomContainer ) {
	    							QueryAtomContainer testQ = (QueryAtomContainer) test;
	    							mcsGhostSMARTSSet.add( smaH.toSmarts( testQ ) );
	    						} else {
	    							try {
										mcsGhostSMARTSSet.add( sGenerator.create( test ) );
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
    	
    	
    	
    	private ExtendedIsomorphism comparator;
    	private MCSMethods mcsMapper;
    	private boolean useWeights;
    	private boolean ghostSubstructures;
    	//private boolean detailedGhostInfo;
    	private int ghostRadius = 4;

    	
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
    	private int uniqueSubstructureCount;
    	private int mcsGhostCount;
    	public List<String> mcsGhostSMARTSSet; 
    	public int modProdNodeCount;
    	public long modProdConstructionTime;
    	public long mcsExecTime;
    	public double modProdEdgeDensity;
    }
    
    
    

}

