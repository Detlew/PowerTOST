# find root of a function to the nearest integer in multiple of step
# code based on ssanv::uniroot.integer, author Michael P. Fay 
uniroot.step <- function (f, interval, lower = min(interval), upper = max(interval), 
                          step=2, step.power = 6, step.up = TRUE, 
                          pos.side = FALSE, maxiter = 100, ...) 
{
  iter <- 0
  if (!is.numeric(lower) || !is.numeric(upper) || lower >= upper) 
    stop("lower < upper  is not fulfilled")
  if (lower == -Inf && step.up == TRUE) 
    stop("lower cannot be -Inf when step.up=TRUE")
  if (upper == Inf && step.up == FALSE) 
    stop("upper cannot be Inf when step.up=FALSE")
  # with step==1 uniroot.integer is obtained
  if (step<1)
    stop("step has to be >=1")
  step <- ceiling(step) # paranoia
  step.power.stop <- 0
  if (step==1) {
    step <- 2
    step.power.stop <- -1
  }  
  
  if (step.up) {
    f.old <- f(lower, ...)
    iter <- iter + 1
    sign <- 1
    xold <- lower
  }
  else {
    f.old <- f(upper, ...)
    iter <- iter + 1
    sign <- -1
    xold <- upper
  }
  ever.switched <- FALSE
  tried.extreme <- FALSE
  while (step.power > step.power.stop) {
    if (f.old == 0) 
      (break)()
    if (iter >= maxiter) 
      stop("reached maxiter without a solution")
    xnew <- xold + sign * step^step.power
    if ((step.up & xnew < upper) || (!step.up & xnew > lower)) {
      f.new <- f(xnew, ...)
      iter <- iter + 1
    }
    else {
      xnew <- xold
      f.new <- f.old
      step.power <- step.power - 1
      if (tried.extreme == FALSE) {
        if (step.up) {
          f.extreme <- f(upper, ...)
          iter <- iter + 1
          x.extreme <- upper
        }
        else {
          f.extreme <- f(lower, ...)
          iter <- iter + 1
          x.extreme <- lower
        }
        tried.extreme <- TRUE
        xswitch <- x.extreme
        f.switch <- f.extreme
        if (f.extreme == 0) {
          xold <- x.extreme
          f.old <- f.extreme
          (break)()
        }
        if (f.old * f.extreme >= 0) {
          stop("f() at extremes not of opposite sign")
        }
      }
    }
    if (f.old * f.new < 0) {
      sign <- sign * (-1)
      ever.switched <- TRUE
      xswitch <- xold
      f.switch <- f.old
    }
    if (ever.switched) {
      step.power <- step.power - 1
      if (step.power == step.power.stop) {
        (break)()
      }
    }
    xold <- xnew
    f.old <- f.new
  }
  if (f.old == 0) {
    root <- xold
    f.root <- f.old
  }
  else if (f.new == 0) {
    root <- xnew
    f.root <- f.new
  }
  else if (f.switch == 0) {
    root <- xswitch
    f.root <- f.switch
  }
  else if (pos.side) {
    # in contrast to uniroot.integer we need '>=' if pos.side=TRUE
    root <- ifelse(f.new >= 0, xnew, xswitch)
    f.root <- ifelse(f.new >= 0, f.new, f.switch)
  }
  else {
    root <- ifelse(f.new < 0, xnew, xswitch)
    f.root <- ifelse(f.new < 0, f.new, f.switch)
  }
  list(iter = iter, f.root = f.root, root = root)
}