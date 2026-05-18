package abyss.lphybeast.tobeast.generator;

import abyss.InformedDirichlet;
import abyss.distributions.InformedDirichletPrior;
import beast.base.core.BEASTInterface;
import beast.base.inference.Distribution;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

/**
 * @author Jasmine Saghafifar
 */
public class InformedDirichletToBEAST implements GeneratorToBEAST<InformedDirichlet, Distribution> {
    @Override
    public Distribution generatorToBEAST(InformedDirichlet generator, BEASTInterface value, BEASTContext context) {
        InformedDirichletPrior beastInformedDirichlet = new InformedDirichletPrior();
        RealVector<PositiveReal> alpha = (RealVector<PositiveReal>) context.getAsRealVector(generator.getConc());
        RealScalar<PositiveReal> scale = (RealScalar<PositiveReal>) context.getAsRealVector(generator.getScale());

        beastInformedDirichlet.setInputValue("alpha", alpha);
        beastInformedDirichlet.setInputValue("scale", scale);
        beastInformedDirichlet.initAndValidate();

        return beastInformedDirichlet;
    }

    @Override
    public Class<InformedDirichlet> getGeneratorClass() {
        return InformedDirichlet.class;
    }

    @Override
    public Class<Distribution> getBEASTClass() {
        return Distribution.class;
    }
}
