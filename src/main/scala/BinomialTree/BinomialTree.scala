package BinomialTree

import MathHelpers.Grid

import scala.math.pow
import scala.util.Random

class BinomialTree(var S0: Double, var r: Double, var sigma: Double, var T: Double, var nbTimes: Int, var nbSimus: Int) {
  private val _dt = T / nbTimes
  private val _times = new Grid(.0, T, nbTimes)
  private val _u = Math.exp(+ sigma * Math.sqrt(_dt))
  private val _d = 1.0 / _u
  private val _pU = (Math.exp(r * _dt) - Math.exp(- sigma * Math.sqrt(_dt))) /
    (Math.exp(+ sigma * Math.sqrt(_dt)) - Math.exp(- sigma * Math.sqrt(_dt)))
  private val _pD = 1.0 - _pU
  private var _grid : Array[Array[Double]] = Array.ofDim[Array[Double]](nbTimes + 1)
  (0 to nbTimes)
    .foreach(iTime => _grid(iTime) = Array.ofDim[Double](iTime + 1))

  private var _paths = Array.ofDim[Double](nbSimus, nbTimes + 1)
  private var _pathIndices = Array.ofDim[Int](nbSimus, nbTimes + 1)

  def GetNbTimes(): Int = {
    return nbTimes
  }

  def Times(): Grid = _times

  def GetGrid(): Array[Array[Double]] = {
    return _grid
  }

  def GetPathsIndices(): Array[Array[Int]] = {
    return _pathIndices
  }

  def S(iT: Int) : Array[Double] = _grid(iT)

  def B(iT: Int): Double = Math.exp(iT * _dt * r)

  def PU() : Double = _pU

  def PD() : Double = _pD

  def RollOut() {
    for(iTime <- 0 to nbTimes) {
      for(jTime <- 0 to iTime) {
        _grid(iTime)(jTime) = S0 * pow(_u, iTime - jTime) * pow(_d, jTime)
      }
    }
  }

  def Simulate(): Unit = {
    RollOut()

    var rand = new Random()

    for(iSimu <- 0 to nbSimus - 1){
      _paths(iSimu)(0) = S0
      _pathIndices(iSimu)(0) = 0

      for(jTime <- 1 to nbTimes){
        var k = 0

        if(rand.nextDouble() <= _pD)
          k = 1;

        _pathIndices(iSimu)(jTime) = _pathIndices(iSimu)(jTime - 1) + k
        _paths(iSimu)(jTime) = _grid(jTime)(_pathIndices(iSimu)(jTime))
      }
    }
  }
}
