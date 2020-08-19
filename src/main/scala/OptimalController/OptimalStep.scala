package OptimalController

class OptimalStep(val t: Double, val Q: Double, val hedge: Double, val nu: Double,
                  val nuIndex : Int, val mu: Double, val S: Double, val B: Double,
                  val nuNext: Int) {
  def GetNuIndex() : Int = nuIndex
  def GetQ() : Double = Q
  def GetHedge() : Double = hedge
  def GetNuNext() : Int = nuNext
  def GetMu() : Double = mu
}

