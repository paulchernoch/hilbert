use criterion::{criterion_group, criterion_main, Criterion};
use std::time::{Duration};
use hilbert::integration::stopping_criteria::StoppingCriteria;
use hilbert::integration::stopping_criteria::StoppingReason;

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

fn estimate_error() {
    let permitted_error = 0.01;
    let mut stopper = StoppingCriteria::new(
        Duration::new(20,0),   // elapsed -> make it unlikely to stop by this criteria
        100,                   // minimum_iterations
        700,                   // permitted_iterations -> make it likely to stop by this criteria (make it larger once we speed things up)
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
    // println!("\nReason:\n    {:?}\n", reason.unwrap());
}


fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("estimate_error", |b| b.iter(|| estimate_error()));
}

criterion_group!(stopping, criterion_benchmark);
criterion_main!(stopping);
