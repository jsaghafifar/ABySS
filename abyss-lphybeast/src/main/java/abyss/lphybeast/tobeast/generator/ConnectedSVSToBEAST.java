package abyss.lphybeast.tobeast.generator;

import abyss.ConnectedSVS;
import abyss.distributions.SVSPrior;
import beast.base.core.BEASTInterface;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

public class ConnectedSVSToBEAST implements GeneratorToBEAST<ConnectedSVS, SVSPrior> {
    @Override
    public SVSPrior generatorToBEAST(ConnectedSVS generator, BEASTInterface value, BEASTContext context) {

        SVSPrior svsPrior = new SVSPrior();
        svsPrior.setInputValue("scale", context.getBEASTObject(generator.getScale()));
        svsPrior.setInputValue("shape", context.getBEASTObject(generator.getShape()));
        svsPrior.setInputValue("rates", context.getBEASTObject(generator.getRates()));
        svsPrior.setInputValue("indicators", value);
        svsPrior.setInputValue("symmetric", generator.getSymmetric().value());
        svsPrior.initAndValidate();
        return svsPrior;
    }

    @Override
    public Class<ConnectedSVS> getGeneratorClass() {
        return ConnectedSVS.class;
    }

    @Override
    public Class<SVSPrior> getBEASTClass() {
        return SVSPrior.class;
    }
}
