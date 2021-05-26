package org.femosaa.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Iterator;
import java.util.Set;

import org.femosaa.core.EAConfigure;
import org.femosaa.core.SASSolution;

import jmetal.core.Solution;
import jmetal.core.SolutionSet;

public class Logger {
	public static String prefix = "/Users/" + System.getProperty("user.name") + "/research/monitor/ws-soa/sas/";
	// This attribute is only used for testing
	public static int max_number_of_eval_to_have_only_seed = 0;

	public static synchronized void logSolutionSet(SolutionSet pareto_front, String name) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = pareto_front.iterator();
			while (itr.hasNext()) {
				Solution s = (Solution) itr.next();
				for (int i = 0; i < s.numberOfObjectives(); i++) {
					if (i == 0) {
						// s.setObjective(0,1.0/(s.getObjective(i) * -1));
					}

					data += s.getObjective(i) + (i == s.numberOfObjectives() - 1 ? "" : ",");
				}
				data += "\n";
			}

			bw.write(data);
			bw.write("------------------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionSetFitness(SolutionSet pareto_front, String name) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = pareto_front.iterator();
			while (itr.hasNext()) {
				Solution s = (Solution) itr.next();
				data += s.getFitness();
				data += "\n";
			}

			bw.write(data);
			bw.write("------------------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logFirstTod(double gen, String name) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = String.valueOf(gen);

			bw.write(data);
			bw.write("------------------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionSetWithGeneration(SolutionSet pareto_front, String name, int gen) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = pareto_front.iterator();
			while (itr.hasNext()) {
				Solution s = (Solution) itr.next();
				for (int i = 0; i < s.numberOfObjectives(); i++) {
					data += s.getObjective(i) + (i == s.numberOfObjectives() - 1 ? "" : ",");
				}
				data += "\n";
			}

			bw.write(data);
			bw.write("------------:" + gen + ":------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logNumberOfNoDominated(SolutionSet pareto_front, String name, int gen) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = String.valueOf(pareto_front.size());
			/*
			 * Iterator itr = pareto_front.iterator(); while(itr.hasNext()) { Solution s =
			 * (Solution)itr.next(); for(int i = 0; i < s.numberOfObjectives(); i++) { data
			 * += s.getObjective(i) + (i == s.numberOfObjectives() - 1? "" : ","); } data +=
			 * "\n"; }
			 */

			bw.write(data + "\n");
			bw.write("------------:" + gen + ":------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logFinalEvaluation(String name, int gen) {
		if (EAConfigure.getInstance().measurement == -1) {
			return;
		}

		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = gen + "\n";

			bw.write(data);
			// bw.write("------------:" + gen + ":------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionSetFitnessWithGeneration(SolutionSet pareto_front, String name,
			int gen) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = pareto_front.iterator();
			while (itr.hasNext()) {
				Solution s = (Solution) itr.next();
				data += s.getFitness();
				data += "\n";
			}

			bw.write(data);
			bw.write("------------:" + gen + ":------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionSetWithGenerationAndFuzzyValue(SolutionSet new_pareto_front,
			SolutionSet old_pareto_front, String name, int gen) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator new_itr = new_pareto_front.iterator();
			Iterator old_itr = old_pareto_front.iterator();
			while (new_itr.hasNext()) {
				Solution s1 = (Solution) new_itr.next();
				Solution s2 = (Solution) old_itr.next();
				boolean isRecorded = true;
				for (int i = 0; i < s1.numberOfObjectives(); i++) {

					if (s2.getObjective(i) == Double.MAX_VALUE || s2.getObjective(i) == (Double.MAX_VALUE / 100)
							|| s1.getObjective(i) == Double.MAX_VALUE
							|| s1.getObjective(i) == (Double.MAX_VALUE / 100)) {
						isRecorded = false;
						break;
					}
				}

				if (!isRecorded) {
					continue;
				}

				for (int i = 0; i < s1.numberOfObjectives(); i++) {
					// original:fuzzy
					data += s2.getObjective(i) + ":" + s1.getObjective(i)
							+ (i == s1.numberOfObjectives() - 1 ? "" : ",");
				}
				data += "\n";
			}

			bw.write(data);
			bw.write("------------:" + gen + ":------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionSetValues(SolutionSet pareto_front, String name) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = pareto_front.iterator();
			while (itr.hasNext()) {
				Solution s = (Solution) itr.next();
				for (int i = 0; i < s.numberOfVariables(); i++) {
					data += s.getDecisionVariables()[i].getValue() + (i == s.numberOfVariables() - 1 ? "" : ",");
				}
				data += "\n";
			}

			bw.write(data);
			bw.write("------------------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionFull(SolutionSet pareto_front, String name, boolean end) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			if (!end) {
				String data = "";
				Iterator itr = pareto_front.iterator();
				while (itr.hasNext()) {
					Solution s = (Solution) itr.next();
					for (int i = 0; i < s.numberOfVariables(); i++) {
						data += s.getDecisionVariables()[i].getValue() + (i == s.numberOfVariables() - 1 ? "" : ":");
					}
					data += "=";

					for (int i = 0; i < s.numberOfObjectives(); i++) {
						data += s.getObjective(i) + (i == s.numberOfObjectives() - 1 ? "" : ",");
					}

					data += "\n";
				}

				bw.write(data);
			} else {
				bw.write("------------------------\n");

			}

			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionFull(Set<String> set, String name) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = set.iterator();
			while (itr.hasNext()) {
				data += itr.next() + "\n";
			}

			bw.write(data);

			bw.write("------------------------\n");

			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionFull(Solution s, String name, boolean end) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			if (!end) {
				String data = "";

				for (int i = 0; i < s.numberOfVariables(); i++) {
					data += s.getDecisionVariables()[i].getValue() + (i == s.numberOfVariables() - 1 ? "" : ":");
				}
				data += "=";

				for (int i = 0; i < s.numberOfObjectives(); i++) {
					data += s.getObjective(i) + (i == s.numberOfObjectives() - 1 ? "" : ",");
				}

				data += "\n";

				bw.write(data);
			} else {
				bw.write("------------------------\n");
			}

			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logSolutionSetValuesWithGen(SolutionSet pareto_front, String name, int gen) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = pareto_front.iterator();
			while (itr.hasNext()) {
				Solution s = (Solution) itr.next();
				for (int i = 0; i < s.numberOfVariables(); i++) {
					data += s.getDecisionVariables()[i].getValue() + (i == s.numberOfVariables() - 1 ? "" : ",");
				}
				data += "\n";
			}

			bw.write(data);
			bw.write("------------:" + gen + ":------------\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void logPercentageOfMarkedSolution(SolutionSet pareto_front, String name) {
		File file = null;
		if (!(file = new File(prefix)).exists()) {
			file.mkdirs();
		}

		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(prefix + name, true));

			String data = "";
			Iterator itr = pareto_front.iterator();
			int no = 0;
			while (itr.hasNext()) {
				Solution s = (Solution) itr.next();
				if (((SASSolution) s).isFromInValid) {
					no++;
				}
			}

			bw.write(no + ":" + no / pareto_front.size() + "\n");
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static synchronized void printMarkedSolution(SolutionSet pareto_front, int eval) {
		Iterator itr = pareto_front.iterator();
		int no = 0;
		while (itr.hasNext()) {
			Solution s = (Solution) itr.next();
			if (((SASSolution) s).isFromInValid) {
				no++;
			}
		}
		System.out.print("from seed: " + no + "\n");
		if (no == pareto_front.size() && max_number_of_eval_to_have_only_seed == 0) {
			max_number_of_eval_to_have_only_seed = eval;
		}
	}
}
