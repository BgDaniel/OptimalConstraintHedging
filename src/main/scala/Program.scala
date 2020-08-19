import BinomialTree.BinomialTree
import OptimalController.OptimalController
import Options.EuropeanCall

object Program {
  def main(args: Array[String]): Unit ={
    ExecuteHedge.Execute();
  }
}

object ExecuteHedge {
  def Execute(): Unit ={
    val S0 = 1.0
    val r = .02
    val sigma = .35
    val T = 1.0
    val nbTimes = 20
    val nbSimus = 1000

    val binomialTree = new BinomialTree(S0, r, sigma, T, nbTimes, nbSimus)
    binomialTree.Simulate()

    val pathIndices = binomialTree.GetPathsIndices()

    val nuSMax = 1.0
    val nuSMin = .0
    val nbNuS = 10
    val VMax = 2.0
    val VMin = -.0
    val nbV = 400

    val penalty = (x: Double, y: Double) => 10.0 * (x - y) * (x - y)

    val transaction = (q: Double) => .03 * Math.abs(q)

    val N = 1.0
    val K = 1.0
    val call = new EuropeanCall(N, K)

    val optimalController = new OptimalController(binomialTree, nuSMin, nuSMax, nbNuS,
      VMin, VMax, nbV, penalty, transaction, call)
    val step0  = optimalController.Control()

    val firstPaths = pathIndices(0)
    var hedge = optimalController.RollOut(firstPaths, step0)
  }

}
