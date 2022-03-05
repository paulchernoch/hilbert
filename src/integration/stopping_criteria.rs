use quantogram::{Quantogram, QuantogramBuilder};
use std::time::{Duration, Instant};

/// Status returned after each iteration of the integration process.
/// These values support the conclusion made in the StoppingReason. 
#[derive(Clone,Copy,Debug)]
pub struct IterationStatus {
    /// Latest estimate of the expected value of the integration.
    /// During the first several score of iterations, this is not reliable.
    pub expected_value: f64, 

    /// Estimated error of the result.
    /// During the first several score of iterations, this is not reliable.
    pub estimated_error: f64, 

    /// expected_value corrected for estimated error.
    pub corrected_value:f64,

    /// Estimated absolute relative error of the result.
    pub absolute_relative_error: f64, 

    /// The deviation is the coefficient of quartile deviation, which 
    /// is (Q3-Q1)/(Q3+Q1), where Q1 and Q3 are the first and third quartile 
    /// values of the estimates.
    pub deviation: f64,

    /// Number of samples received so far (including NAN and infinity values).
    pub iterations: u64,

    /// Time elapsed since the first sample was received.
    pub elapsed: Duration
}

/// Gives the reason that the integration should be halted (if it should), 
/// along with the values used to justify the decision, such as the 
/// latest estimates of the expected value and error. 
#[derive(Clone,Copy,Debug)]
pub enum StoppingReason {
    /// The execution time limit was reached. 
    /// The result may not be accurate. 
    TimeLimit(IterationStatus),

    /// The maximum number of samples was reached. 
    /// The result may not be accurate.
    SampleLimit(IterationStatus),

    /// The estimated absolute relative error is small enough.
    LowEstimatedError(IterationStatus),

    /// The Coefficient of Quartile Deviation is small enough,
    /// meaning most estimates are confined to a narrow range.
    SmallDeviation(IterationStatus),

    /// No reason to stop yet.
    Continue(IterationStatus)
}

impl StoppingReason {
    pub fn should_stop(&self) -> bool {
        match self {
            StoppingReason::Continue(_) => false,
            _ => true
        }
    }
    pub fn get_status(&self) -> IterationStatus {
        match self {
            StoppingReason::TimeLimit(s) => s.clone(),
            StoppingReason::SampleLimit(s) => s.clone(),
            StoppingReason::LowEstimatedError(s) => s.clone(),
            StoppingReason::SmallDeviation(s) => s.clone(),
            StoppingReason::Continue(s) => s.clone()
        }
    }
}

/// Decides when integration should stop based on multiple criteria.
/// 
/// 1. Out of time.
/// 2. Exceeds maximum iteration budget.
/// 3. Estimated error is low enough. 
/// 4. Series of estimated solutions has a strong enough peak around the mode.
/// 
/// Every iteration supplies three values:
/// 
///   - latest estimate for expected value
///   - change to previous estimate
///   - latest sample value
/// 
/// If a sample is NaN, the estimate and change will be disregarded, but the value will count
/// toward the iteration budget.
/// 
/// The typical use of this class is to call `evaluate` in a loop and quit when 
/// a value other than StoppingReason::Continue(status) is return.
pub struct StoppingCriteria {
    /// Time budget for the integration.
    permitted_elapsed_time: Duration,

    /// Barring a timeout, perform at least this many 
    /// iterations before checking the other criteria.
    minimum_iterations: u64,

    /// Iteration budget for the integration. 
    /// NAN and infinity values count toward this value.
    permitted_iterations: u64,

    /// Permitted absolute relative error. 
    /// Stop if the estimated error drops below this.
    /// NOTE: Quantogram defines a best accuracy of 0.0347%. 
    /// This cannot do better than that. 
    permitted_error: f64,

    /// Permitted value of the coefficient of quartile deviation 
    /// (a value between zero and one) for the estimates.
    /// A small value indicates a sharp peak of estimates, meaning that 
    /// most estimates fall within a narrow band, hence have settled upon a value.
    /// A good value for this will be similar in magnitude to permitted_error.
    /// If estimates come in asymptotic from above or below, it may
    /// be advisiable to make this half the value of permitted_error,
    /// as the estimate is taken from the mode, and in those cases the mode will
    /// be off center relative to the peak.
    permitted_coeff_of_quartile_dev: f64,

    /// Histogram of all estimates of the expected value.
    /// If the sampled value is NAN, the corresponding estimate is not added to the histogram.
    /// The final guess will be the half-sample mode of the estimates.
    estimates: Quantogram,

    /// Histogram of all changes to the estimate.
    /// If the sampled value is NAN, the corresponding change is not added to the histogram.
    deltas: Quantogram,

    /// Count of how many samples have been taken that are not NAN or infinity.
    finite_sample_count: u64,

    /// Count of how many samples have been taken including NAN or infinity.
    all_sample_count: u64,

    /// Time When the first sample was recorded.
    start_time: Option<Instant>,

    last_error_check: usize

}

impl StoppingCriteria {
    /// Create a StoppingCriteria to monitor progress towards a solution. 
    ///   
    /// # Arguments
    /// 
    /// * `permitted_elapsed_time` - Direct algorithm to stop if this much execution time is exceeded.
    /// * `minimum_iterations` - Prevent algorithm frmo stopping until at least this many iterations occur. 
    /// * `permitted_iterations` - Direct algorithm to stop when this many iterations have occurred.
    /// * `permitted_error` - Direct algorithm to stop when an estimate of the absolute relative error drops below this value.
    /// * `permitted_coeff_of_quartile_dev` - Direct algorithm to stop when the series of estimates cluster together tightly enough that
    /// the coefficient of quartile deviation drops below this value.
    pub fn new(
        permitted_elapsed_time: Duration,
        minimum_iterations: u64,
        permitted_iterations: u64,
        permitted_error: f64,
        permitted_coeff_of_quartile_dev: f64
    ) -> Self {
        StoppingCriteria {
            permitted_elapsed_time: permitted_elapsed_time,
            minimum_iterations: minimum_iterations,
            permitted_iterations: permitted_iterations,
            permitted_error: permitted_error,
            permitted_coeff_of_quartile_dev: permitted_coeff_of_quartile_dev,
            // The Quantogram abs rel error must be smaller than the desired error for the integration. 
            // Making it 1/10th slowed things down a lot due to hsm calcs. 1/2 works well.
            estimates: QuantogramBuilder::new().with_error(permitted_error / 2.0).build(),
            deltas: QuantogramBuilder::new().with_error(permitted_error / 2.0).build(),
            finite_sample_count: 0,
            all_sample_count: 0,
            start_time: None,
            last_error_check: 0
        }
    }

    /// Incorporate a new sample, updated estimate and its change since the previous estimate.
    /// Return an indicator of whether it is time to stop, why we should stop and whether we have
    /// a good estimate or not.
    /// 
    /// # Arguments
    /// 
    /// * `sample` - An unweighted sample of the integrand. 
    /// * `estimate` - Latest estimate of the expected value of the integration.
    /// * `delta` - Change to the estimate since the previous iteration.
    pub fn evaluate(&mut self, sample: f64, estimate: f64, delta:f64) -> StoppingReason {
        // Evaluation of the first sample starts the clock.
        if self.start_time.is_none() {
            self.start_time = Some(Instant::now());
        }

        self.all_sample_count += 1;
        if sample.is_finite() {
            self.finite_sample_count += 1;
            self.estimates.add(estimate);
            self.deltas.add(delta);
        }

        // Make sure some error checks are performed at nice round numbers of iterations.
        let do_checks_at: Vec<usize> = vec![
            1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 
            15000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 
            200000, 500000, 1000000, 2000000, 5000000, 10000000
        ];

        let check_error = if self.all_sample_count < 100 {
            // Always check the estimated error early on.
            self.last_error_check = self.all_sample_count as usize;
            true 
        }
        else if self.all_sample_count as f64 >= self.last_error_check as f64 * 1.095 || do_checks_at.contains(&(self.all_sample_count as usize)) {
            // Spread out error checks logarithmically (they are costly).
            self.last_error_check = self.all_sample_count as usize;
            true 
        }
        else { 
            false
        };

        let (error, abs_rel_error, exp_value) = if check_error {self.err_absrelerr_result() } else { (f64::NAN, f64::NAN, f64::NAN) };

        let mut status = IterationStatus {
            expected_value: exp_value, 
            estimated_error: error, 
            corrected_value: exp_value - error,
            absolute_relative_error: abs_rel_error,
            deviation: self.deviation(), 
            iterations: self.all_sample_count, 
            elapsed: self.elapsed_time()
        };

        let mut reason = self.choose_reason(status);

        if reason.should_stop() && !check_error {
            // If there is a timeout or we exceed iterations, then 
            // make sure the caller gets a full assessment, so we need to check the error after all. 
            let (error, abs_rel_error, exp_value) = self.err_absrelerr_result();
            status.expected_value = exp_value;
            status.estimated_error = error;
            status.corrected_value = exp_value - error;
            status.absolute_relative_error = abs_rel_error;
            reason = self.choose_reason(status);
        }

        reason
    }

    /// Estimate the error as -2N x mode(deltas). 
    /// This is based on the following assumptions: 
    ///   - Error drops proportional to 1/√N for Monte Carlo integration (random sampling)
    ///   - Errors for Pseudo-random sampling are usually better and no worse than for random sampling
    ///   - We are taking the limit as N -> ∞
    ///   - The mode of the series of deltas is a good smoothing function
    ///   - half-sample mode is a good estimate of the mode for histogrammed values
    /// 
    /// Returns NAN if no finite values have been seen yet.
    pub fn estimate_error(&self) -> f64 {
        // The negative sign is because if the deltas are positive, it means we think we have to add 
        // more to get to the true value, so our error is negative. 
        match self.deltas.hsm() {
            Some(mode) => -2.0 * self.finite_sample_count as f64 * mode,
            None => f64::NAN
        }
    }

    /// Estimate the expected value as mode(estimates). 
    /// 
    /// This is based on the following assumptions: 
    ///   - As estimates converge on a solution, they form a peak in the distribution
    ///   - The half-sample mode is a good measure of that peak
    /// 
    /// Returns NAN if no finite values have been seen yet.
    pub fn estimate_result(&self) -> f64 {
        match self.estimates.hsm() {
            Some(mode) => mode,
            None => f64::NAN
        }
    }

    /// Compute the estimated error, absolute relative error and expected value
    /// all at once, to avoid repeating expensive calculations.
    fn err_absrelerr_result(&self) -> (f64, f64, f64) {
        let error = self.estimate_error();
        let value = self.estimate_result();
        let abs_rel_error = if error.is_finite() && value.is_finite() && value != 0.0 {
            (error/value).abs()
        }
        else {
            f64::NAN
        };
        (error, abs_rel_error, value)
    }

    /// Estimate the absolute relative error.
    pub fn absolute_relative_error(&self) -> f64 {
        let (_error, abs_rel_error, _result) = self.err_absrelerr_result();
        abs_rel_error
    }

    /// Get the coefficient of quartile deviation for the estimate, which in terms
    /// of the 1st and 3rd quartiles is (Q3-Q1)/(Q3+Q1).
    /// 
    /// The value is between zero and one. 
    /// A low value indicates a sharp peak and likely convergence.
    /// For example, if this deviation is 0.02 (2%), if the mode is at the center of
    /// the peak, then that means half of the estimates fall within ± 1% of the estimate
    /// for the expected value. If estimates are asymptotic from below or above,
    /// the mode will be off-center, and only a quarter of the estimates are within 1%. 
    pub fn deviation(&self) -> f64 {
        match self.estimates.coeff_of_quartile_dev() {
            Some(coeff) => coeff,
            None => f64::NAN
        }
    }

    /// Decide if and why integration should stop. 
    /// If multiple reasons for stopping are detected, choose the most desirable,
    /// in this order:
    ///   - LowEstimatedError
    ///   - SmallDeviation
    ///   - SampleLimit
    ///   - TimeLimit
    fn choose_reason(&self, status: IterationStatus) -> StoppingReason {
        if status.iterations < self.minimum_iterations {
            StoppingReason::Continue(status)
        }
        else if status.absolute_relative_error <= self.permitted_error {
            StoppingReason::LowEstimatedError(status)
        }
        else if status.deviation <= self.permitted_coeff_of_quartile_dev {
            StoppingReason::SmallDeviation(status)
        }
        else if status.iterations >= self.permitted_iterations {
            StoppingReason::SampleLimit(status)
        }
        else if status.elapsed >= self.permitted_elapsed_time {
            StoppingReason::TimeLimit(status)
        }
        else {
            StoppingReason::Continue(status)
        }      
    }

    /// Time elapsed since first sample was evaluated.
    pub fn elapsed_time(&self) -> Duration {
        match self.start_time {
            Some(t) => t.elapsed(),
            None => Duration::ZERO
        } 
    }
}



#[cfg(test)]
mod tests {
    use std::time::{Duration};
    use super::StoppingCriteria;
    use super::StoppingReason;

    /// Test that when convergence criteria are evaluated against a simulated
    /// series of successive approximations, it concludes that the error is low enough
    /// after the expected number of iterations and with an error estimate that is in the right ballpark.
    /// 
    /// The series behaves as though the initially very high 100% error drops as `1/√N` according to theory. 
    /// When that error is assumed to be below 1%, we stop. 
    /// 
    /// Since the initial error is 100, after 10,000 iterations, it shoud drop by 1/sqrt(10000) or 1/100. 
    #[test]
    fn estimate_error() {
        let permitted_error = 0.01; // Permit 1% error.
        let mut stopper = StoppingCriteria::new(
            Duration::new(40,0),   // elapsed -> make it unlikely to stop by this criteria
            100,                   // minimum_iterations
            100000,                // permitted_iterations -> make it unlikely to stop by this criteria
            permitted_error,       // permitted_error -> hope that abs rel error is the reason for stopping
            0.0                    // permitted_coeff_of_quartile_dev -> make it impossible by this criteria
        );
        let true_value = 100.0;
        let init_error = 100.0;
        let mut reason: Option<StoppingReason> = None;
        let mut iteration: usize = 1;
        let mut prev_estimate = 0.0;
        while reason.is_none() || !reason.unwrap().should_stop() {
            let estimate = estimate_with_error(iteration, true_value, init_error, 0.7);
            let sample_doesnt_matter = 0.0;
            let delta = estimate - prev_estimate;
            reason = Some(stopper.evaluate(sample_doesnt_matter, estimate, delta));
            prev_estimate = estimate;
            iteration += 1;
        }
        println!("\nReason:\n    {:?}\n", reason.unwrap());
        match reason {
            Some(StoppingReason::LowEstimatedError(status)) => {
                let buest_guess = status.corrected_value;
                let actual_abs_rel_error = ((buest_guess - true_value) / true_value).abs();
                assert!(actual_abs_rel_error <= permitted_error);
                // It actually takes about 10000 iterations. 
                // Since we do not check the error after each iteration, just be in the right ballpark.
                assert!(status.iterations <= 11000);
            },
            _ => assert!(false)
        }
    }

    /// Golden ratio based pseudo-random number sequence.
    ///  ```text
    ///    φ = (1 + √5) / 2
    ///    G(n) = (n · φ) mod 1
    /// ```
    fn golden_ratio_sequence(n: usize) -> f64 {
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        (n as f64 * phi) % 1.0
    }

    /// Simulate a series of estimates of a value where the error is explicitly given 
    /// and its reduction in magnitude is defined to be proportional to `1/√N` where N is the number of iterations.
    /// 
    /// # Arguments
    /// * `iteration` - Place in the sequence, which should be positive (start from one).
    /// * `true_value` - True value of the quantity being estimated.
    /// * `go_high_frequency` - A fraction [0,1] that indicates how often the estimate should be high and how often low. 
    /// If one, always go high. If zero, always go low. 
    /// * returns the next estimate for the value. 
    fn estimate_with_error(iteration: usize, true_value: f64, initial_error: f64, go_high_frequency: f64) -> f64 {
        let pseudo_random = golden_ratio_sequence(iteration);
        let sign: f64 = 
            if go_high_frequency >= 1.0 { 1.0 }
            else if go_high_frequency <= 0.0 { -1.0 }
            else if pseudo_random < go_high_frequency { 1.0 }
            else { -1.0 };
        let current_error = sign * initial_error / (iteration as f64).sqrt();
        true_value + current_error
    }

}
