package abyss;

import cern.colt.bitvector.BitVector;

/**
 * ported from beast1 to BEAST_CLASSIC - author: Marc Suchard
 * @author dkuh004
 *         Date: Sep 18, 2011
 *         Time: 6:14:07 PM
 */
public interface SVS {

    class Utils {

        public static boolean connectedAndWellConditioned(double[] probability) {
            for(double prob : probability) {
                if(prob < tolerance || prob >= 1.0) {
                    return false;
                }
            }
            return true;
        }

        public static void setTolerance(double newTolerance) {
            tolerance = newTolerance;
        }

        public static double getTolerance() {
            return tolerance;
        }

        /* Determines if the graph is strongly connected, such that there exists
        * a directed path from any vertex to any other vertex
        *
        */
        public static boolean isStronglyConnected(Boolean[] indicatorValues, int dim, boolean reversible) {
            BitVector visited = new BitVector(dim);
            boolean connected = true;
            for (int i = 0; i < dim && connected; i++) { //TODO optimise so doesn't have to check connectivity from every node
                visited.clear();
                depthFirstSearch(i, visited, indicatorValues, dim, reversible);
                connected = visited.cardinality() == dim;
            }
            return connected;
        }

        private static boolean hasEdge(int i, int j, Boolean[] indicatorValues,
                                       int dim, boolean reversible) {
            return i != j && indicatorValues[getEntry(i, j, dim, reversible)];
        }

        private static int getEntry(int i, int j, int dim, boolean reversible) {
            if (reversible) {
                if (j < i) {
                    return getEntry(j,i,dim,reversible);
                }
                int entry = i * dim - i * (i + 1) / 2 + j - 1 -i;
                return entry;
            }

            int entry = i * (dim - 1) + j;
            if (j > i)
                entry--;
            return entry;
        }

        private static void depthFirstSearch(int node, BitVector visited, Boolean[] indicatorValues,
                                             int dim, boolean reversible) {
            visited.set(node);
            for (int v = 0; v < dim; v++) {
                if (hasEdge(node, v, indicatorValues, dim, reversible) && !visited.get(v))
                    depthFirstSearch(v, visited, indicatorValues, dim, reversible);
            }
        }

        private static double defaultExpectedMutations = 1.0;
        private static double tolerance = 1E-20;
    }
    
}   // interface BSSVS

