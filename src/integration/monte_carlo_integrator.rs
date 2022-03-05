use rand::prelude::*;
use super::stopping_criteria::{StoppingCriteria,StoppingReason};

/// Set a random point of the appropriate number of dimensions D with the given lower and upper bounds.
pub fn random_point<const D: usize>(p: &mut [f64; D], lower_bounds: [f64; D], upper_bounds: [f64; D], rng: &mut ThreadRng) {
    for dimension in 0..D {
        p[dimension] = rng.gen_range(lower_bounds[dimension]..upper_bounds[dimension]);
    }
}

fn golden_ratio_ranged(n: usize, low: f64, high: f64) -> f64 {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let r = (n as f64 * phi) % 1.0;
    low + (high - low) * r
}

pub fn pseudo_random_point<const D: usize>(p: &mut [f64; D], lower_bounds: [f64; D], upper_bounds: [f64; D], n: &mut usize) {
    for dimension in 0..D {
        p[dimension] = golden_ratio_ranged(*n, lower_bounds[dimension], upper_bounds[dimension]);
        *n += 1;
    }
}

/// Perform Monte-Carlo integration of an integrand function in D dimensions.
pub fn integrate_monte_carlo<F, const D: usize> (integrand: F, lower_bounds: [f64; D], upper_bounds: [f64; D], pseudo_random: bool, stopper: &mut StoppingCriteria) -> StoppingReason
where F: Fn(&[f64]) -> f64 {
    let mut sum: f64 = 0.0;
    let mut iteration: usize = 0;
    let mut area = 1.0;
    let mut point: [f64; D] = [0.0; D];
    let mut rng = rand::thread_rng();
    for (low, high) in lower_bounds.iter().zip(upper_bounds.iter()) {
        area *= high - low;
    }
    let mut reason: Option<StoppingReason> = None;
    let mut prev_estimate = 0.0;
    let mut n = 0_usize;
    while reason.is_none() || !reason.unwrap().should_stop() {
        if pseudo_random {
            pseudo_random_point(&mut point, lower_bounds, upper_bounds, &mut n);
        }
        else {
            random_point(&mut point, lower_bounds, upper_bounds, &mut rng);
        }
        let sample = integrand(&point);
        sum += sample;
        iteration += 1;
        let estimate = area * sum / (iteration as f64);
        let delta = estimate - prev_estimate;
        reason = Some(stopper.evaluate(sample, estimate, delta));
        prev_estimate = estimate;
    }
    // println!("\n\n   *** FINAL ESTIMATE: {:.6}  ***", estimate);
    reason.unwrap()
}



#[cfg(test)]
mod tests {
    use std::time::{Duration};
    use super::StoppingCriteria;
    use super::integrate_monte_carlo;

    #[test]
    fn monte_carlo_2dims() {
        let permitted_error = 0.001; // Permit 0.1% error.
        let mut stopper = StoppingCriteria::new(
            Duration::new(20,0),   // elapsed 
            100,                   // minimum_iterations
            100000,                // permitted_iterations
            permitted_error,       // permitted_error 
            0.0005                 // permitted_coeff_of_quartile_dev
        );
        let f = |p: &[f64]| -> f64 { 
            let x = p[0];
            let y = p[1];
            let z = 4.0 + (10.0 * x * y * (x+y).sin() / (x*x + y*y + 1.0));
            z
        };
        let lower: [f64;2] = [-8.0,-8.0];
        let upper: [f64;2] = [8.0,8.0];
        // In tests, pseudo random is giving smaller error with 1/20th as many iterations 
        // compared to true random.
        let pseudo_random = true;
        let reason = integrate_monte_carlo(f, lower, upper, pseudo_random, &mut stopper);
        let status = reason.get_status();
        let reason_string = format!("\n\n  Reason = {:?}\n\n", reason)
            .replace(",", "\n    ")
            .replace("{ ", "{\n     ")
            .replace("})","\n})");
        println!("{}", reason_string);
        let true_value = 1024.0;
        let true_abs_relative_error = (status.corrected_value - true_value).abs() / true_value;
        assert!(true_abs_relative_error < 0.001);
        // assert!(false);
    }

}