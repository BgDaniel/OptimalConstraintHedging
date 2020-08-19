package Options

class EuropeanCall(val N: Double, val K: Double) extends TerminalOption {
  override def Pay(S: Double): Double = N * Math.max(S - K, .0)
}
