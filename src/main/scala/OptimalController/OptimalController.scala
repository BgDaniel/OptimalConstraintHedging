package OptimalController

import java.util

import BinomialTree.BinomialTree
import MathHelpers.Grid
import Options.TerminalOption

class OptimalController(var binomialTree: BinomialTree, var nuMin: Double, var nuMax: Double, var nbNu: Int,
                        var VMin: Double, var VMax: Double, var nbV: Int,
                        var penalty: (Double, Double) => Double, var transaction: Double => Double,
                        var option: TerminalOption) {

  private val _nu = new Grid(nuMin, nuMax, nbNu)
  private val _V = new Grid(VMin, VMax, nbV)
  private val _nbTimes = binomialTree.GetNbTimes()
  private val _grid = binomialTree.GetGrid()
  private val _times = binomialTree.Times()

  private var _optimalSteps : Array[Array[Array[Array[OptimalStep]]]] =
    Array.ofDim[Array[OptimalStep]](_nbTimes + 1, nbV + 1, nbNu + 1)

  for(t <- 0 to _nbTimes) {
    for(v <- 0 to nbV) {
      for (n <- 0 to nbNu) {
        _optimalSteps(t)(v)(n) = Array.ofDim[OptimalStep](_grid(t).length)
      }
    }
  }

  for(iV <- 0 to nbV) {
    for (jNu <- 0 to nbNu) {
      for(kS <- 0 to _nbTimes){
        val t = _times.Val(_nbTimes)
        val hedge = _V.Val(iV)
        val nu = _nu.Val(jNu)
        val S = binomialTree.S(_nbTimes)(kS)
        val B = binomialTree.B(_nbTimes)
        val mu = (hedge - nu * S) / B
        val Q = penalty(hedge, option.Pay(S))
        _optimalSteps(_nbTimes)(iV)(jNu)(kS) = new OptimalStep(t, Q, hedge, nu, jNu, mu, S, B, -1)
      }
    }
  }

  def DetermineQ(iT: Int, jV: Int, kNu: Int, lS: Int): (Double, Int) ={
    var _QNew = Double.MaxValue
    var QNew = .0
    var nuNext = 0

    for(lNu <- 0 to nbNu){
      var deltaNu = _nu.Val(lNu) - _nu.Val(kNu)
      var transactionCosts = transaction(deltaNu)
      var expectedQ = binomialTree.PU() * _optimalSteps(iT + 1)(jV)(lNu)(lS).GetQ() +
                      binomialTree.PD() * _optimalSteps(iT + 1)(jV)(lNu)(lS + 1).GetQ()
      var QNew = transactionCosts + expectedQ

      if(QNew < _QNew){
        _QNew = QNew
        nuNext = lNu
      }
    }

    return (_QNew, nuNext)
  }

  def RollOut(pathIndices : Array[Int], step0: OptimalStep) : (Array[Double], Array[Double], Array[Double]) ={
    val Q : Array[Double] =  Array.ofDim[Double](_nbTimes)
    val hedge : Array[Double] =  Array.ofDim[Double](_nbTimes)
    val S : Array[Double] =  Array.ofDim[Double](_nbTimes)

    Q(0) = step0.GetQ()
    hedge(0) = step0.GetHedge()
    S(0) = binomialTree.S(0)(pathIndices(0))

    var nuIndexNext = step0.GetNuIndex()
    var optimalStep = step0

    for(iT <- 1 to _nbTimes - 1) {
      var nuNext = _nu.Val(nuIndexNext)
      val muNext = (hedge(iT - 1) - nuNext * S(iT - 1)) / binomialTree.B(iT - 1)
      S(iT) = binomialTree.S(iT)(pathIndices(iT))
      hedge(iT) = nuNext * S(iT) + muNext * binomialTree.B(iT)

      val iV = _V.Idx(hedge(iT))
      optimalStep = _optimalSteps(iT)(iV)(nuIndexNext)(pathIndices(iT))
      Q(iT) = optimalStep.GetQ()

      nuIndexNext = optimalStep.GetNuNext()
    }

    return (Q, hedge, S)
  }

  def Control(): OptimalStep = {
    for(iT <- _nbTimes - 1 to 0 by -1) {
      for(jV <- 0 to nbV) {
        for (kNu <- 0 to nbNu) {
          for(lS <- 0 to iT){
            val t = _times.Val(iT)
            val hedge = _V.Val(jV)
            val nu = _nu.Val(kNu)
            val S = binomialTree.S(iT)(lS)
            val B = binomialTree.B(iT)
            val mu = (hedge - nu * S) / B

            val v = DetermineQ(iT, jV, kNu, lS)
            _optimalSteps(iT)(jV)(kNu)(lS) = new OptimalStep(t, v._1, hedge, nu, kNu, mu, S, B, v._2)
          }
        }
      }
    }

    var _QMin = Double.MaxValue
    var step0 = new OptimalStep(0, 0, 0, 0, 0, 0, 0, 0, -1)

    for(jV <- 0 to nbV) {
      for (kNu <- 0 to nbNu) {
        val Q = _optimalSteps(0)(jV)(kNu)(0).GetQ()

        if (Q < _QMin)
          {
            _QMin = Q
            step0 = _optimalSteps(0)(jV)(kNu)(0)
          }
      }
    }

    return step0
  }
}
