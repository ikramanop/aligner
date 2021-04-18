use std::str;

#[derive(Debug, Clone, Copy)]
pub enum Direction {
    Up,
    Left,
    Diagonal,
    Beginning,
}

#[derive(Debug, Clone, Copy)]
pub enum Nucleotide {
    A = 64,
    G = 71,
    C = 67,
    T = 84,
}

#[derive(Debug, Clone, Copy)]
pub enum Protein {
    A = 0,
    R = 1,
    N = 2,
    D = 3,
    C = 4,
    Q = 5,
    E = 6,
    G = 7,
    H = 8,
    I = 9,
    L = 10,
    K = 11,
    M = 12,
    F = 13,
    P = 14,
    S = 15,
    T = 16,
    W = 17,
    Y = 18,
    V = 19,
    B = 20,
    J = 21,
    Z = 22,
    X = 23,
    Any = 24,
}

impl Protein {
    pub fn match_with_u8(symbol: u8) -> Protein {
        let result = match symbol {
            65 => Protein::A,
            82 => Protein::R,
            78 => Protein::N,
            68 => Protein::D,
            67 => Protein::C,
            81 => Protein::Q,
            69 => Protein::E,
            71 => Protein::G,
            72 => Protein::H,
            73 => Protein::I,
            76 => Protein::L,
            75 => Protein::K,
            77 => Protein::M,
            70 => Protein::F,
            80 => Protein::P,
            83 => Protein::S,
            84 => Protein::T,
            87 => Protein::W,
            89 => Protein::Y,
            86 => Protein::V,
            66 => Protein::B,
            74 => Protein::J,
            90 => Protein::Z,
            88 => Protein::X,
            42 => Protein::Any,
            s => panic!("No such protein {:?}", str::from_utf8(&[s])),
        };

        result
    }

    pub fn u8_vec_to_protein_vec(sequence: &[u8]) -> Vec<Protein> {
        let mut result = Vec::<Protein>::new();

        for elem in sequence.iter() {
            result.push(Protein::match_with_u8(*elem));
        }

        result
    }
}
