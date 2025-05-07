package abyss.spi;

import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;
import lphy.base.spi.LPhyBaseImpl;
import abyss.*;

import java.util.Arrays;
import java.util.List;

/**
 *
 * The provider of SPI which is an implementation of a service.
 * It requires a public no-args constructor.
 * @author Walter Xie
 */
public class ABySSImpl extends LPhyBaseImpl {
    /**
     * Required by ServiceLoader.
     */
    public ABySSImpl() {
        //TODO print package or classes info here?
    }

    @Override
    public List<Class<? extends GenerativeDistribution>> declareDistributions() {
        return Arrays.asList( );
    }

    @Override
    public List<Class<? extends BasicFunction>> declareFunctions() {
        return Arrays.asList(NonReversible.class, NQPFAM.class);
    }

    public String getExtensionName() {
        return "Advanced BaYesian Site and Substitution models";
    }
}
