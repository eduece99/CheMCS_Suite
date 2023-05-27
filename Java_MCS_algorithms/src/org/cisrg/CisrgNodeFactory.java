package org.cisrg;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Cisrg" Node.
 * 
 *
 * @author EdD
 */
public class CisrgNodeFactory 
        extends NodeFactory<CisrgNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public CisrgNodeModel createNodeModel() {
        return new CisrgNodeModel();
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
    public NodeView<CisrgNodeModel> createNodeView(final int viewIndex,
            final CisrgNodeModel nodeModel) {
        return new CisrgNodeView(nodeModel);
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
        return new CisrgNodeDialog();
    }

}

