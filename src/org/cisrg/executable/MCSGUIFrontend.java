package org.cisrg.executable;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.Border;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JTabbedPane;
import javax.swing.JFileChooser;


import org.cisrg.knime.GraphSimilarityNodeModel.AggregationMethod;
import org.cisrg.knime.GraphSimilarityNodeModel.ReturnType;
import org.cisrg.knime.GraphSimilarityNodeModel.SimilarityType;
import org.cisrg.mapping.ConvenienceTools;
import org.cisrg.mapping.ExtendedAlgorithm;
import org.cisrg.mapping.SimilarityComparator;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import picocli.CommandLine.Option;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;



public class MCSGUIFrontend extends Frame {
	
	
	private void addMoleculesToPanel( List<IAtomContainer> mols ) {
		
		simHub.calculateSimilarity(mols.get(0), mols.get(1));
		
		Image bi = null;
		DepictionGenerator dptgen = new DepictionGenerator();
		
		GridBagConstraints gbc = new GridBagConstraints();
	    gbc.gridx = 0;
        gbc.gridy = 0;
		
		for( int mi = 0; mi < mols.size(); mi++ ) {
			
			IAtomContainer mol = mols.get(mi);
			
		    try {
		    	  dptgen = dptgen.withSize(200, 250)              // px (raster) or mm (vector)
				      .withMolTitle()
				      .withTitleColor(Color.DARK_GRAY); // annotations are red by default

				  
				  if( mi == 0 ) {
					  dptgen = dptgen.withHighlight( simHub.bondMaps.get(0).keySet() , Color.RED).withOuterGlowHighlight(2.0) ;
				  } else {
					  dptgen = dptgen.withHighlight( simHub.bondMaps.get(0).values() , Color.RED).withOuterGlowHighlight(2.0) ;
				  }
				  
				  bi = dptgen.depict(mol).toImg();
							
			} catch (CDKException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		    
		    
	        
	        molPanel.add( new JLabel( new ImageIcon(bi) ), gbc );
	        gbc.gridx++;
		}
		
		molPanel.repaint();
		molPanel.validate();
	}
	
	
	class FileButtonListener implements ActionListener {

		JButton buttonComponent;

		public FileButtonListener( JButton bc ) {
			super();
			
			this.buttonComponent = bc;
		}
		
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			
			currentButton = buttonComponent;
			fc.showOpenDialog(buttonComponent);
			
		}
		
	}
	
	class LoadMoleculesButtonListener implements ActionListener {

		public LoadMoleculesButtonListener() {
			super();
		}
		
		@Override
		public void actionPerformed(ActionEvent ev) {
			// TODO Auto-generated method stub
			
			// get molecules
	        ArrayList<IAtomContainer> compounds = null;
			try {
				refMols = ConvenienceTools.getQueryMolecules( new File( refFileComponent.getText() ), null, true ); 
				dbMols = ConvenienceTools.getQueryMolecules( new File( dbFileComponent.getText() ), null, true );
				//System.out.println( "compound file size - " + compounds.size() );
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				
				JOptionPane.showMessageDialog(null,"Input molecule files do not exist.  Check your file paths","Error",1);
				
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				
				JOptionPane.showMessageDialog(null,"Cannot read input molecules.  Only .smi and .sdf Files are supported","Error",1);
				
				e.printStackTrace();
			}
			
			
			// remove all molecules
			molPanel.removeAll();
			
			List<IAtomContainer> mols = new ArrayList<IAtomContainer>(2);
			mols.add( refMols.get(0) );
			mols.add( dbMols.get(0) );
			
			addMoleculesToPanel( mols );
			
		}
		
	}
	
	class FileChooserListener implements ActionListener {
		
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			
			if( currentButton == refFileButtonComponent ) {
				refFileComponent.setText( fc.getSelectedFile().getAbsolutePath() );
			} else if( currentButton == dbFileButtonComponent ) {
				dbFileComponent.setText( fc.getSelectedFile().getAbsolutePath() );
			} else {
				System.out.println("No file path saved." + currentButton );
				
				
			}
		}
		
	}

	
	JComboBox mappingAlgorithmComponent;
	JComboBox aggregationMethodComponent;
	JComboBox returnTypeComponent;
	JComboBox similarityTypeComponent;
	
	JTextField topoDistanceLimitComponent;
	JTextField timeLimitComponent;
	JTextField refFileComponent;
	JTextField dbFileComponent;
	JButton refFileButtonComponent;
	JButton dbFileButtonComponent;
	JButton loadMoleculesButtonComponent;
	
	JButton currentButton;
	
	JPanel molPanel;
	
	JCheckBox outputMCS;
	JCheckBox raymondHeuristicsComponent;
	JCheckBox ringHeuristicsComponent;
	JCheckBox detailedGhostInfoComponent;
	
	List<IAtomContainer> refMols, dbMols;
	
	final JFileChooser fc;
	
	SimilarityComparator simHub = null;
	private String algorithmName = "Depolli_dMCES";
    private int topologicalDistanceLimit = 4;
    private boolean ringHeuristics = false;
    private boolean raymondHeuristics = false;
    private int timeLimit = 10000;
	private boolean ghostSubstructures = false;
	private boolean bondWeightFlag = false;
	
	
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
	
	
	/** Enum for the similarity measure types. */
	public enum SimilarityType { 
		Tanimoto, Tversky, Number
	}

	public MCSGUIFrontend() {
		
		/*
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
       */
		
		
		ExtendedAlgorithm algorithm = ExtendedAlgorithm.valueOf(algorithmName);
		
		simHub = new SimilarityComparator(
				null, bondWeightFlag, ghostSubstructures, algorithm, 
				raymondHeuristics, ringHeuristics, topologicalDistanceLimit, timeLimit, false
		);
		
		refFileComponent = new JTextField("input file path", 200);  
		refFileButtonComponent = new JButton("Browse");
		refFileButtonComponent.addActionListener( new FileButtonListener(refFileButtonComponent) );
		
		dbFileComponent = new JTextField("input file path", 200); 
		dbFileButtonComponent = new JButton("Browse");
		dbFileButtonComponent.addActionListener( new FileButtonListener(dbFileButtonComponent) );
		
		
		loadMoleculesButtonComponent = new JButton("Load Molecules");
		loadMoleculesButtonComponent.addActionListener( new LoadMoleculesButtonListener() );
		
		
		fc = new JFileChooser();
		fc.addActionListener( new FileChooserListener( ) );
		fc.addActionListener( new FileChooserListener( ) );
       
        mappingAlgorithmComponent = new JComboBox( );
        mappingAlgorithmComponent.setModel( new DefaultComboBoxModel( MappingAlgorithmNames  ) );
        mappingAlgorithmComponent.setSelectedIndex(1);
        
        aggregationMethodComponent = new JComboBox();
        aggregationMethodComponent.setModel( new DefaultComboBoxModel( AggregationMethod.values() ) );
        
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
        
        molPanel = new JPanel() ;
        molPanel.setLayout(new GridBagLayout() );
        molPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        
        JTabbedPane tp = new JTabbedPane();  
        //tp.setLayout( new GridBagLayout()  );
        tp.setBounds(50,50,200,200);  
        
        JPanel ioPanel = new JPanel();  
        ioPanel.setLayout(new GridBagLayout() );
        
        JPanel mcsSettingsPanel = new JPanel();  
        mcsSettingsPanel.setLayout(new GridBagLayout() );
        
        
        
        tp.add("IO", ioPanel);
        tp.add("MCS", mcsSettingsPanel);
        
        
        
        //addTab("Settings", mainPanel);
        
        gbc.fill = GridBagConstraints.BOTH;
        gbc.gridx = 0;
        gbc.gridy = 0;
        ioPanel.add(new JLabel("reference CDK Column: "), gbc);
        gbc.gridy++;
        
        ioPanel.add(new JLabel("database CDK Column: "), gbc);
        gbc.gridy++;
        
        gbc.gridx = 1;
        gbc.gridy = 0;
        
        gbc.fill = GridBagConstraints.BOTH ;
        gbc.weightx = 0.8;
        ioPanel.add( refFileComponent, gbc);
        gbc.gridy++;
        
        ioPanel.add( dbFileComponent, gbc);
        gbc.gridy++;
        
        gbc.weightx = 0.3;
        
        gbc.gridx = 2;
        gbc.gridy = 0;
        
        ioPanel.add( refFileButtonComponent, gbc);
        gbc.gridy++;
        
        ioPanel.add( dbFileButtonComponent, gbc);
        gbc.gridy++;
        
        
        gbc.gridx = 0;
        //gbc.gridy = 0;
        ioPanel.add( loadMoleculesButtonComponent, gbc);
        gbc.gridy++;
        
        
        mcsSettingsPanel.add(new JLabel("search algorithm: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("aggregation method: "), gbc);
        gbc.gridy++;
        
        //mainPanel.add(new JLabel("return method: "), gbc);
        //gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("similarity equation: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("Modular Product heuristics: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("Modular Product Ring heuristics: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("topological distance limit: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("expansion time limit: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("output MCS: "), gbc);
        gbc.gridy++;
        
        mcsSettingsPanel.add(new JLabel("output ghost SMARTS: "), gbc);
        gbc.gridy++;
        
        
        
        gbc.gridx = 1;
        gbc.gridy = 0;
        //mainPanel.add( refMolColumnComponent, gbc );
        gbc.gridy++;
        
        //mainPanel.add( databaseMolColumnComponent, gbc );
        gbc.gridy++;
             
        mcsSettingsPanel.add( mappingAlgorithmComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( aggregationMethodComponent, gbc );
        gbc.gridy++;
        
        //mainPanel.add( returnTypeComponent, gbc );
        //gbc.gridy++;
        
        mcsSettingsPanel.add( similarityTypeComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( raymondHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( ringHeuristicsComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( topoDistanceLimitComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( timeLimitComponent, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( outputMCS, gbc );
        gbc.gridy++;
        
        mcsSettingsPanel.add( detailedGhostInfoComponent, gbc );
        gbc.gridy++;
		
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent windowEvent){
               System.exit(0);
            }        
         });    
        
        
        IChemObjectBuilder bldr   = SilentChemObjectBuilder.getInstance();
        SmilesParser       smipar = new SmilesParser(bldr);
        
        
        IAtomContainer mol = null;
        BufferedImage bi = null;
		try {
			mol = smipar.parseSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine");
		} catch (InvalidSmilesException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        mol.setProperty(CDKConstants.TITLE, "caffeine"); // title already set from input!

        DepictionGenerator dptgen = new DepictionGenerator();
        try {
			 bi = dptgen.withSize(200, 250)              // px (raster) or mm (vector)
			      .withMolTitle()
			      .withTitleColor(Color.DARK_GRAY) // annotations are red by default
			      .depict(mol).toImg();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        gbc.gridx = 0;
        gbc.gridy = 0;
        molPanel.setBackground(Color.red);
        molPanel.add( new JLabel( new ImageIcon(bi) ), gbc );
        gbc.gridy++;
              
        
        setSize(600,400);  
        setLayout( new GridBagLayout() );
        
        
        // add panels to Frame
        
        gbc.gridx = 0;
        gbc.gridy = 0;
        add( molPanel , gbc );
        
        gbc.gridx = 0;
        gbc.gridy = 1;
        add( mainPanel , gbc );
        
        
        gbc.gridx = 0;
        gbc.gridy = 2;
        add( tp , gbc );
        
        setTitle("CheMCS GUI");
		setVisible(true);  
	}
	
	
	
	public static void main(String[] args) {
		
		try {
            // Set System L&F
			//UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
			UIManager.setLookAndFeel("com.sun.java.swing.plaf.gtk.GTKLookAndFeel");
		} 
	    catch (UnsupportedLookAndFeelException e) {
	       e.printStackTrace();
	       
	    }
	    catch (ClassNotFoundException e) {
	    	e.printStackTrace();
	    }
	    catch (InstantiationException e) {
	    	e.printStackTrace();
	    }
	    catch (IllegalAccessException e) {
	    	e.printStackTrace();
	    }
		
		MCSGUIFrontend f = new MCSGUIFrontend();    

	}
	
	

}
