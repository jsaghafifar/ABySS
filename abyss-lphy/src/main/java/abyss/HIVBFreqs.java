package abyss;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

//TODO model citation
public class HIVBFreqs<U> extends DeterministicFunction<U[]> {

    public HIVBFreqs() {}

    @GeneratorInfo(name = "hivbFreqs", description = "")
    public Value<U[]> apply() {
        HIVB hivb = new HIVB(null);
        U[] array = (U[]) hivb.getFreq();
        return new Value<>( array, this);
    }
}