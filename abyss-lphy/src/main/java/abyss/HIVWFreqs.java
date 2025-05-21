package abyss;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

//TODO model citation
public class HIVWFreqs<U> extends DeterministicFunction<U[]> {

    public HIVWFreqs() {}

    @GeneratorInfo(name = "hivwFreqs", description = "")
    public Value<U[]> apply() {
        HIVW hivw = new HIVW(null);
        U[] array = (U[]) hivw.getFreq();
        return new Value<>( array, this);
    }
}