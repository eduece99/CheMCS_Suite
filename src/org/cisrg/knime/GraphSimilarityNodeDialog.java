package org.cisrg.knime;

import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.Arrays;
import java.util.Vector;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.Border;

import org.cisrg.knime.GraphSimilarityNodeModel.AggregationMethod;
import org.cisrg.knime.GraphSimilarityNodeModel.ReturnType;
import org.cisrg.knime.GraphSimilarityNodeModel.SimilarityType;
import org.cisrg.mapping.ExtendedAlgorithm;
import org.knime.chem.types.SmartsValue;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.DataValue;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.util.ColumnSelectionPanel;
import org.knime.core.node.util.DataValueColumnFilter;
import org.openscience.cdk.knime.type.CDKValue;
import org.openscience.cdk.knime.commons.CDKNodeUtils;

/**
 * <code>NodeDialog</code> for the "GraphSimilarity" Node.
 * Graph-based similarity searching with support for hyperstructures
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author CISRG (Edmund Duesbury)
 */
public class GraphSimilarityNodeDialog extends NodeDialogPane {

    /**
     * New pane for configuring GraphSimilarity node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
    protected GraphSimilarityNodeDialog() {
    super();
     
    	Class<? extends DataValue>[] accepted = new Class[ CDKNodeUtils.ACCEPTED_VALUE_CLASSES.length+1 ];
    	System.arraycopy(CDKNodeUtils.ACCEPTED_VALUE_CLASSES, 0, accepted, 0, CDKNodeUtils.ACCEPTED_VALUE_CLASSES.length);
    	accepted[ accepted.length-1 ] = SmartsValue.class;
 
 
    	
        refMolColumnComponent = new ColumnSelectionPanel( 
    			(Border) null, 
    			//new DataValueColumnFilter( StringValue.class ) ,
    			new DataValueColumnFilter( CDKValue.class, SmartsValue.class ) ,
    			//new DataValueColumnFilter( accepted ) ,
    			false  // dummy column if true
    	);
        
        databaseMolColumnComponent = new ColumnSelectionPanel( 
    			(Border) null, 
    			new DataValueColumnFilter( CDKValue.class ) ,
    			//new DataValueColumnFilter( CDKNodeUtils.ACCEPTED_VALUE_CLASSES ) ,
    			false  // dummy column if true
    	);
       
        mappingAlgorithmComponent = new JComboBox( );
        mappingAlgorithmComponent.setModel( new DefaultComboBoxModel( GraphSimilarityNodeModel.MappingAlgorithmNames  ) );
        mappingAlgorithmComponent.setSelectedIndex(1);
        
        aggregationMethodComponent = new JComboBox();
        aggregationMethodComponent.setModel( new DefaultComboBoxModel( AggregationMethod.values() ) );
        
        returnTypeComponent = new JComboBox();
        returnTypeComponent.setModel( new DefaultComboBoxModel( ReturnType.values() ) );
        returnTypeComponent.setSelectedIndex(0);
        
        similarityTypeComponent = new JComboBox();
        similarityTypeComponent.setModel( new DefaultComboBoxModel( SimilarityType.values() ) );
        similarityTypeComponent.setSelectedIndex(0);
        
        outputMCS = new JCheckBox();
        outputMCS.setSelected(false);
        
        topoDistanceLimitComponent = new JTextField(2);  // 2 columns
        topoDistanceLimitComponent.setText("-1");
        
        timeLimitComponent = new JTextField(6);  // 6 columns
        timeLimitComponent.setText("10000");
        
        raymondHeuristicsComponent = new JCheckBox();
        raymondHeuristicsComponent.setSelected(false);
        
        ringHeuristicsComponent = new JCheckBox();
        ringHeuristicsComponent.setSelected(false);
        
        detailedGhostInfoComponent = new JCheckBox(); 
        detailedGhostInfoComponent.setSelected(false);
        		
        
        GridBagConstraints gbc = new GridBagConstraints();
        
        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new GridBagLayout() );
        mainPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        
        
        addTab("Settings", mainPanel);
        
        gbc.fill = GridBagConstraints.BOTH;
        gbc.gridx = 0;
        gbc.gridy = 0;
        mainPanel.add(new JLabel("reference CDK Column: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("database CDK Column: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("search algorithm: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("aggregation method: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("return method: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("similarity equation: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("Modular Product heuristics: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("Modular Product Ring heuristics: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("topological distance limit: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("expansion time limit: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("output MCS: "), gbc);
        gbc.gridy++;
        
        mainPanel.add(new JLabel("output ghost SMARTS: "), gbc);
        gbc.gridy++;
        
        
        
        gbc.gridx = 1;
        gbc.gridy = 0;
        mainPanel.add( refMolColumnComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( databaseMolColumnComponent, gbc );
        gbc.gridy++;
             
        mainPanel.add( mappingAlgorithmComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( aggregationMethodComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( returnTypeComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( similarityTypeComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( raymondHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( ringHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( topoDistanceLimitComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( timeLimitComponent, gbc );
        gbc.gridy++;
        
        mainPanel.add( outputMCS, gbc );
        gbc.gridy++;
        
        mainPanel.add( detailedGhostInfoComponent, gbc );
        gbc.gridy++;
    }

    
	@Override
	protected void saveSettingsTo(NodeSettingsWO settings)
			throws InvalidSettingsException {
		
		settings.addString(GraphSimilarityNodeModel.CFGKEY_REFMOLCOLUMN, refMolColumnComponent.getSelectedColumn() );
		settings.addString(GraphSimilarityNodeModel.CFGKEY_DATABASEMOLCOLUMN, databaseMolColumnComponent.getSelectedColumn() );
		
		
		settings.addString(GraphSimilarityNodeModel.CFGKEY_ALGORITHM, mappingAlgorithmComponent.getSelectedItem().toString() );
		settings.addString(GraphSimilarityNodeModel.CFGKEY_AGGREGATION, aggregationMethodComponent.getSelectedItem().toString() );
		settings.addString(GraphSimilarityNodeModel.CFGKEY_RETURNTYPE, returnTypeComponent.getSelectedItem().toString() );
		settings.addString(GraphSimilarityNodeModel.CFGKEY_SIMILARITYTYPE, similarityTypeComponent.getSelectedItem().toString() );
		
		settings.addBoolean(GraphSimilarityNodeModel.CFGKEY_OUTPUTMCS, outputMCS.isSelected() );
		settings.addBoolean(GraphSimilarityNodeModel.CFGKEY_RAYMONDHEURISTICS, raymondHeuristicsComponent.isSelected() );
		settings.addBoolean(GraphSimilarityNodeModel.CFGKEY_RINGHEURISTICS, ringHeuristicsComponent.isSelected() );
		settings.addBoolean(GraphSimilarityNodeModel.CFGKEY_MOREGHOSTINFO, detailedGhostInfoComponent.isSelected() );
		
		// XXX error could occur here (non-number entry)
		int tdl = Integer.parseInt( topoDistanceLimitComponent.getText() );
		settings.addInt(GraphSimilarityNodeModel.CFGKEY_TOPODISTANCECONSTRAINT, tdl);
		
		// XXX error could occur here (non-number entry)
		String tls = timeLimitComponent.getText();
		if( tls.length() > 6 )
			tls = tls.substring(0,6);
		int tl = Integer.parseInt( tls );
		settings.addInt(GraphSimilarityNodeModel.CFGKEY_EXPTIMEOUT, tl);
	}
	
	
	@Override
	protected void loadSettingsFrom(NodeSettingsRO settings,
			DataTableSpec[] specs) throws NotConfigurableException {
		
		try {

			
			databaseMolColumnComponent.update(specs[0], settings.getString(GraphSimilarityNodeModel.CFGKEY_DATABASEMOLCOLUMN) );
			refMolColumnComponent.update(specs[1], settings.getString(GraphSimilarityNodeModel.CFGKEY_REFMOLCOLUMN) );
			//System.out.println( settings.getDataType( GraphSimilarityNodeModel.CFGKEY_MOLCOLUMN ) );
			//molColumnComponent.setSelectedColumn( settings.getString(HyperstructureConstructorNodeModel.CFGKEY_MOLCOLUMN) );
			
			mappingAlgorithmComponent.setSelectedItem( 
					ExtendedAlgorithm.valueOf(
							settings.getString(GraphSimilarityNodeModel.CFGKEY_ALGORITHM)
					)
			);
			aggregationMethodComponent.setSelectedItem( 
					GraphSimilarityNodeModel.AggregationMethod.valueOf( 
							settings.getString(GraphSimilarityNodeModel.CFGKEY_AGGREGATION) 
					) 
			);
			returnTypeComponent.setSelectedItem( 
					GraphSimilarityNodeModel.ReturnType.valueOf( 
							settings.getString(GraphSimilarityNodeModel.CFGKEY_RETURNTYPE ) 
					) 
			);
			similarityTypeComponent.setSelectedItem( 
					GraphSimilarityNodeModel.SimilarityType.valueOf( 
							settings.getString(GraphSimilarityNodeModel.CFGKEY_SIMILARITYTYPE ) 
					) 
			);
			
			outputMCS.setSelected( settings.getBoolean(GraphSimilarityNodeModel.CFGKEY_OUTPUTMCS ) );
			raymondHeuristicsComponent.setSelected( settings.getBoolean(GraphSimilarityNodeModel.CFGKEY_RAYMONDHEURISTICS) );
			ringHeuristicsComponent.setSelected( settings.getBoolean(GraphSimilarityNodeModel.CFGKEY_RINGHEURISTICS) );
			detailedGhostInfoComponent.setSelected( settings.getBoolean(GraphSimilarityNodeModel.CFGKEY_MOREGHOSTINFO ) );
			
			topoDistanceLimitComponent.setText( settings.getInt( GraphSimilarityNodeModel.CFGKEY_TOPODISTANCECONSTRAINT ) + "" );
			timeLimitComponent.setText( settings.getInt( GraphSimilarityNodeModel.CFGKEY_EXPTIMEOUT ) + "" );
			
			/*
			if ( GraphSimilarityNodeModel.aggregationMethod == AggregationMethod.Minimum) {
				m_minimum.setSelected(true);
			} else if ( GraphSimilarityNodeModel.aggregationMethod == AggregationMethod.Maximum) {
				m_maximum.setSelected(true);
			} else if ( GraphSimilarityNodeModel.aggregationMethod == AggregationMethod.Average) {
				m_average.setSelected(true);
			} else if ( GraphSimilarityNodeModel.aggregationMethod == AggregationMethod.Matrix) {
				m_matrix.setSelected(true);
			}

			if (m_settings.returnType().equals(ReturnType.String)) {
				returnString.setSelected(true);
			} else if (m_settings.returnType().equals(ReturnType.Collection)) {
				returnCollection.setSelected(true);
			}

			if (m_settings.aggregationMethod() == AggregationMethod.Matrix) {
				returnString.setEnabled(false);
				returnCollection.setEnabled(false);
			} else {
				returnString.setEnabled(true);
				returnCollection.setEnabled(true);
			}
			*/
			
		} catch (InvalidSettingsException e) {
			e.printStackTrace();
			System.out.println("Settings missing in load method!");
		}
	}
			
	/*
	private final ColumnSelectionComboxBox molColumnComponent = new ColumnSelectionComboxBox( 
			(Border) null, 
			DataValue.class
	);
	*/
	ColumnSelectionPanel refMolColumnComponent;
	ColumnSelectionPanel databaseMolColumnComponent;
	
	JComboBox mappingAlgorithmComponent;
	JComboBox aggregationMethodComponent;
	JComboBox returnTypeComponent;
	JComboBox similarityTypeComponent;
	
	JTextField topoDistanceLimitComponent;
	JTextField timeLimitComponent;
	
	JCheckBox outputMCS;
	JCheckBox raymondHeuristicsComponent;
	JCheckBox ringHeuristicsComponent;
	JCheckBox detailedGhostInfoComponent;
}

