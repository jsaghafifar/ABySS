package abyss.lphybeast.tobeast.generator;

import abyss.InformedDirichlet;
import abyss.distributions.InformedDirichletPrior;
import beast.base.core.BEASTInterface;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

/**
 * @author Jasmine Saghafifar
 */
public class InformedDirichletToBEAST implements GeneratorToBEAST<InformedDirichlet, Prior> {
    @Override
    public Prior generatorToBEAST(InformedDirichlet generator, BEASTInterface value, BEASTContext context) {
        InformedDirichletPrior beastInformedDirichlet = new InformedDirichletPrior();
        Value<Number[]> alpha = generator.getConc();
        Value<Number> scale = generator.getScale();

        beastInformedDirichlet.setInputValue("alpha", context.getBEASTObject(alpha));
        beastInformedDirichlet.setInputValue("scale", context.getBEASTObject(scale));
        beastInformedDirichlet.initAndValidate();

        return BEASTContext.createPrior(beastInformedDirichlet, (RealParameter) value);
    }

    @Override
    public Class<InformedDirichlet> getGeneratorClass() {
        return InformedDirichlet.class;
    }

    @Override
    public Class<Prior> getBEASTClass() {
        return Prior.class;
    }
}
