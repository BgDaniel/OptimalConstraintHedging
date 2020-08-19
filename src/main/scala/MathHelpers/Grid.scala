package MathHelpers

class Grid(min: Double, max: Double, nbSteps: Int) {
  private val _dx = (max - min) / nbSteps
  private val _x: Array[Double] = Array.tabulate(nbSteps + 1)(i => min + i * _dx)

  def Idx(v: Double) : Int = {
    if (v < min)
      throw new Exception(String.format("Value %s smaller than lower grid boundary %d!", v, min))

    if (v > max)
      throw new Exception(String.format("Value %s greater than upper grid boundary %d!", v, max))

    val ell = Math.floor((v - min) / _dx).toInt

    if (v - ell <= .5)
      return ell
    else
      return ell + 1
  }

  def Val(ell: Int) : Double = _x(ell)
}
