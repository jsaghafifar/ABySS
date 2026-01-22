package abyss;

import lphy.base.evolution.substitutionmodel.RateMatrix;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

//TODO citation here
public class FLU extends RateMatrix {

    public FLU(@ParameterInfo(name = RateMatrix.meanRateParamName, description = "the mean rate of the process. default 1.0", optional=true) Value<Number> meanRate) {
        super(meanRate);
    }

    @GeneratorInfo(name = "flu", verbClause = "is", narrativeName = "Empirical influenza amino acid substitution model",
            category = GeneratorCategory.RATE_MATRIX,
            description = "A symmetric instantaneous rate matrix for influenza amino acid substitution.")
    public Value<Double[][]> apply() {
        Double[][] Q = getQ();
        return new DoubleArray2DValue(Q, this);
    }

    protected Double[][] getQ() {
        Double[] r = getRates();
        Double[] f = getFreq();
        int numStates = f.length;
        Double[][] Q = new Double[numStates][numStates];
        int x = 0;
        for (int i = 0; i < f.length; i++) {
            for (int j = i + 1; j < f.length; j++) {
                Q[i][j] = r[x] * f[j];
                Q[j][i] = r[x] * f[i];
                x++;
            }
            Q[i][i] = 0.0;
            for (int j = 0; j < f.length; j++) {
                if (i==j) continue;
                Q[i][i] -= Q[i][j];
            }
        }
        normalize(f, Q, totalRateDefault1());
        return Q;
    }

    public Double[] getRates() {
        Double[] rates;
        rates = new Double[]{0.026447095,
                             0.584852306,0.014100000,
                             1.484234503,0.000000000,5.370511279,
                             0.080490909,0.104053666,0.000001060,0.001003501,
                             1.132313122,0.116941459,1.934832784,1.593098825,0.001236645,
                             0.214757862,0.021800000,0.887570549,0.256491863,0.119028506,0.058774527,
                             0.149926734,0.001112158,0.014085917,0.014200000,1.463357278,0.000016300,0.243190142,
                             0.474333610,0.000003830,0.290042980,3.881488809,0.320000000,0.264148929,0.347302791,0.227707997,
                             0.023116952,0.005613627,0.005730682,0.016499536,2.986800036,0.006516229,0.321611694,3.512072282,0.129223639,
                             0.058745423,0.111457310,0.041762964,0.313974351,0.279910509,0.001500467,0.001273509,9.017954203,1.331291619,6.746936485,
                             0.053366579,0.000013000,7.737392871,0.061652192,0.000836000,0.322524648,1.387096032,0.218571975,2.646847965,0.000836000,0.005251688,
                             0.659311478,0.000000000,0.188539456,0.319558828,0.007132430,0.038631761,0.924466914,0.080543327,0.195750632,0.634308521,0.056900000,0.036400000,
                             0.353753982,0.002547334,0.145469388,1.195629122,0.032680657,0.108051341,5.330313412,0.028839950,2.559587177,1.020366955,0.190259181,0.530642655,0.712769599,
                             0.138658765,0.167207008,0.006771843,0.124897617,0.016100000,1.190624465,1.879569938,0.246117172,15.30009662,0.296045557,0.890162346,0.161000889,0.154027180,3.292716942,
                             3.011344519,0.336263345,0.338372183,0.307140298,0.996685670,1.585646577,0.580704250,0.290381075,0.283807672,0.570766693,0.007026588,3.881310531,2.087385344,0.487822499,0.950138410,
                             5.418298175,0.011975266,0.135481233,0.280124895,0.000134906,0.018808030,0.368713573,2.904052286,1.526964200,0.044926357,2.031511321,2.140332316,0.542251094,0.602340963,0.183076905,2.206859934,
                             3.532005270,0.054904564,0.297123975,0.285047948,0.592587985,0.337229619,0.098631355,14.39405219,0.073127930,0.890598579,4.904842235,0.010257517,0.058971975,0.406697814,0.103964386,0.088256423,0.654109108,
                             0.196000000,0.094106680,0.000014900,0.155245492,0.814753094,0.196486447,0.022400000,0.032132150,0.000049800,0.431277663,0.070460039,0.000536000,0.000431021,0.044000000,1.369429408,0.099835753,0.207066206,0.256900461,
                             0.018289288,0.601692431,0.525398543,0.104092870,5.393924245,0.074814997,6.448954446,0.273934263,0.012416222,0.340058468,0.874272175,0.373101927,0.000182000,0.072205935,0.099855497,0.392552240,0.124898020,0.167581647,0.42775543};
        double sum = 0;
        for (int i = 0; i < rates.length; i++) {
            sum += rates[i];
        }
        for (int i = 0; i < rates.length; i++) {
            rates[i] /= sum;
        }
        return rates;
    }

    public Double[] getFreq() {
        Double[] freq = new Double[20]; // ORDER: ACDEFGHIKLMNPQRSTVWY
        freq[0] = 0.04707195;
        freq[1] = 0.02502197;
        freq[2] = 0.04785995;
        freq[3] = 0.05458695;
        freq[4] = 0.03049597;
        freq[5] = 0.07637292;
        freq[6] = 0.01996398;
        freq[7] = 0.06713393;
        freq[8] = 0.05678494;
        freq[9] = 0.07149793;
        freq[10]= 0.01815098;
        freq[11]= 0.07421393;
        freq[12]= 0.05065595;
        freq[13]= 0.03330397;
        freq[14]= 0.05090995;
        freq[15]= 0.08840891;
        freq[16]= 0.07433893;
        freq[17]= 0.06322894;
        freq[18]= 0.01852398;
        freq[19]= 0.03147397;
        double sum = 0;
        for (int i = 0; i < freq.length; i++) {
            sum += freq[i];
        }
        for (int i = 0; i < freq.length; i++) {
            freq[i] /= sum;
        }
        return freq;
    }

}