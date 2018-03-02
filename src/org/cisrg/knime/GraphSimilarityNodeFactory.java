package org.cisrg.knime;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GraphSimilarity" Node.
 * Graph-based similarity searching with support for hyperstructures
 *
 * @author CISRG (Edmund Duesbury)
 */
public class GraphSimilarityNodeFactory 
        extends NodeFactory<GraphSimilarityNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GraphSimilarityNodeModel createNodeModel() {
        return new GraphSimilarityNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<GraphSimilarityNodeModel> createNodeView(final int viewIndex,
            final GraphSimilarityNodeModel nodeModel) {
        return new GraphSimilarityNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new GraphSimilarityNodeDialog();
    }

}

