use crate::engine::task::Task;
use crate::filter;
use aligner_core::alignment::PWMAlignment;

#[test]
fn filter_test() {
    let alignment = PWMAlignment::empty();

    let tasks = vec![
        Task {
            alignment: alignment.clone(),
            z: 12.240966,
            left_coord: 300,
            right_coord: 630,
        },
        Task {
            alignment: alignment.clone(),
            z: 12.378159,
            left_coord: 360,
            right_coord: 690,
        },
        Task {
            alignment: alignment.clone(),
            z: 11.762683,
            left_coord: 1080,
            right_coord: 1410,
        },
        Task {
            alignment: alignment.clone(),
            z: 10.471823,
            left_coord: 1740,
            right_coord: 2070,
        },
        Task {
            alignment: alignment.clone(),
            z: 11.392030,
            left_coord: 1860,
            right_coord: 2190,
        },
    ];

    let expected_result = vec![
        Task {
            alignment: alignment.clone(),
            z: 12.378159,
            left_coord: 360,
            right_coord: 690,
        },
        Task {
            alignment: alignment.clone(),
            z: 11.762683,
            left_coord: 1080,
            right_coord: 1410,
        },
        Task {
            alignment,
            z: 11.392030,
            left_coord: 1860,
            right_coord: 2190,
        },
    ];

    assert_eq!(filter(tasks).unwrap(), expected_result)
}
