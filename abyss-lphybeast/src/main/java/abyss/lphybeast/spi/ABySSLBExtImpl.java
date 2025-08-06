package abyss.lphybeast.spi;

import abyss.*;
import beast.base.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.core.model.Generator;
import lphy.core.vectorization.operation.Slice;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;
import abyss.lphybeast.tobeast.generator.*;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The "Container" provider class of SPI
 * which include a list of {@link ValueToBEAST},
 * {@link GeneratorToBEAST}, and {@link DataType}
 * to extend.
 * @author Walter Xie
 */
public class ABySSLBExtImpl implements LPhyBEASTExt {
    // the first matching converter is used.
    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return Arrays.asList( );
    }

    // the first matching converter is used.
    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList( NonReversibleToBEAST.class, NQPFAMToBEAST.class,
                InformedDirichletToBEAST.class, ConnectedSVSToBEAST.class,
                MixedAlignmentToBEAST.class);
    }

    // LPhy SequenceType => BEAST DataType
    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        return new ConcurrentHashMap<>();
    }

    //*** these below are extra from Exclusion, only implemented in extensions ***//

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return List.of(NQPFAMRates.class, NQPFAMFreqs.class,
                FLURates.class, FLUFreqs.class,
                HIVBRates.class, HIVBFreqs.class,
                HIVWRates.class, HIVWFreqs.class,
                ComputeEquilibrium.class, Slice.class);
    }

    @Override
    public List<Class> getExcludedValueType() {
        // For a complex logic, or arrays, use isExcludedValue
        return List.of( );
    }
}
