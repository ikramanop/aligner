use crate::{Error, Result};
use ndarray::Array1;
use rand::{thread_rng, Rng};
use std::cmp::PartialEq;
use std::fmt::Debug;
use std::hash::Hash;
use std::str;

#[derive(Debug, Clone, Copy)]
pub enum Direction {
    Top,
    Left,
    Diagonal,
    Beginning,
}

impl Direction {
    pub(crate) fn get_direction(top: f64, left: f64, diagonal: f64) -> (f64, Direction) {
        let max = f64::max(f64::max(top, left), diagonal);

        if (max - top).abs() < f64::EPSILON {
            (max, Direction::Top)
        } else if (max - left).abs() < f64::EPSILON {
            (max, Direction::Left)
        } else {
            (max, Direction::Diagonal)
        }
    }

    pub(crate) fn get_direction_with_beginning(
        top: f64,
        left: f64,
        diagonal: f64,
    ) -> (f64, Direction) {
        let max = f64::max(f64::max(top, left), diagonal);

        if max == 0f64 {
            (max, Direction::Beginning)
        } else if (max - top).abs() < f64::EPSILON {
            (max, Direction::Top)
        } else if (max - left).abs() < f64::EPSILON {
            (max, Direction::Left)
        } else {
            (max, Direction::Diagonal)
        }
    }
}

impl PartialEq for Direction {
    fn eq(&self, other: &Self) -> bool {
        *self as usize == *other as usize
    }
}

#[derive(Debug, Clone, Copy, Eq)]
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
    Blank = 98,
    Pos = 99,
    Any,
}

impl Hash for Protein {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        core::mem::discriminant(self).hash(state);
    }
}

impl PartialEq for Protein {
    fn eq(&self, other: &Self) -> bool {
        *self as usize == *other as usize
    }
}

impl From<Protein> for usize {
    fn from(p: Protein) -> Self {
        p as usize
    }
}

impl From<usize> for Protein {
    fn from(d: usize) -> Self {
        match d {
            0 => Protein::A,
            1 => Protein::R,
            2 => Protein::N,
            3 => Protein::D,
            4 => Protein::C,
            5 => Protein::Q,
            6 => Protein::E,
            7 => Protein::G,
            8 => Protein::H,
            9 => Protein::I,
            10 => Protein::L,
            11 => Protein::K,
            12 => Protein::M,
            13 => Protein::F,
            14 => Protein::P,
            15 => Protein::S,
            16 => Protein::T,
            17 => Protein::W,
            18 => Protein::Y,
            19 => Protein::V,
            20 => Protein::B,
            21 => Protein::J,
            22 => Protein::Z,
            23 => Protein::X,
            98 => Protein::Blank,
            99 => Protein::Pos,
            _ => Protein::Any,
        }
    }
}

#[derive(Debug, Clone, Copy, Eq)]
pub enum DNA {
    A = 0,
    T = 1,
    C = 2,
    G = 3,
    Blank = 98,
    Pos = 99,
    Any,
}

impl From<DNA> for usize {
    fn from(d: DNA) -> Self {
        d as usize
    }
}

impl From<usize> for DNA {
    fn from(d: usize) -> Self {
        match d {
            0 => DNA::A,
            1 => DNA::T,
            2 => DNA::C,
            3 => DNA::G,
            98 => DNA::Blank,
            99 => DNA::Pos,
            _ => DNA::Any,
        }
    }
}

impl Hash for DNA {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        core::mem::discriminant(self).hash(state);
    }
}

impl PartialEq for DNA {
    fn eq(&self, other: &Self) -> bool {
        *self as usize == *other as usize
    }
}

pub trait BioData
where
    Self: Sized,
{
    fn match_with_char(symbol: char) -> Result<Self>;
    fn convert_to_char(elem: &Self) -> Result<char>;
    fn str_to_vec(sequence: &str) -> Result<Vec<Self>>;
    fn vec_to_str(sequence: &[Self]) -> Result<String>;
    fn from_u8_vec(vec: &[u8]) -> Result<Vec<Self>>;
    fn from_u8_vec_with_freqs(vec: &[u8]) -> Result<(Vec<Self>, Array1<f64>)>;
    fn from_u8_vec_with_freqs_and_indices(
        vec: &[u8],
    ) -> Result<(Vec<Self>, Array1<f64>, Vec<Index>)>;
    fn random_seq(length: usize) -> Result<Vec<Self>>;
    fn random_seq_with_freqs(length: usize) -> Result<(Vec<Self>, Array1<f64>)>;
    fn blank() -> Self;
    fn pos() -> Self;
    fn volume() -> usize;
}

impl BioData for Protein {
    fn match_with_char(symbol: char) -> Result<Protein> {
        match symbol {
            'A' => Ok(Protein::A),
            'R' => Ok(Protein::R),
            'N' => Ok(Protein::N),
            'D' => Ok(Protein::D),
            'C' => Ok(Protein::C),
            'Q' => Ok(Protein::Q),
            'E' => Ok(Protein::E),
            'G' => Ok(Protein::G),
            'H' => Ok(Protein::H),
            'I' => Ok(Protein::I),
            'L' => Ok(Protein::L),
            'K' => Ok(Protein::K),
            'M' => Ok(Protein::M),
            'F' => Ok(Protein::F),
            'P' => Ok(Protein::P),
            'S' => Ok(Protein::S),
            'T' => Ok(Protein::T),
            'W' => Ok(Protein::W),
            'Y' => Ok(Protein::Y),
            'V' => Ok(Protein::V),
            'B' => Ok(Protein::B),
            'J' => Ok(Protein::J),
            'Z' => Ok(Protein::Z),
            'X' => Ok(Protein::X),
            '_' => Ok(Protein::Blank),
            '+' => Ok(Protein::Pos),
            _ => Err(Error::CharIsNotMatchable),
        }
    }

    fn convert_to_char(protein: &Protein) -> Result<char> {
        match protein {
            Protein::A => Ok('A'),
            Protein::R => Ok('R'),
            Protein::N => Ok('N'),
            Protein::D => Ok('D'),
            Protein::C => Ok('C'),
            Protein::Q => Ok('Q'),
            Protein::E => Ok('E'),
            Protein::G => Ok('G'),
            Protein::H => Ok('H'),
            Protein::I => Ok('I'),
            Protein::L => Ok('L'),
            Protein::K => Ok('K'),
            Protein::M => Ok('M'),
            Protein::F => Ok('F'),
            Protein::P => Ok('P'),
            Protein::S => Ok('S'),
            Protein::T => Ok('T'),
            Protein::W => Ok('W'),
            Protein::Y => Ok('Y'),
            Protein::V => Ok('V'),
            Protein::B => Ok('B'),
            Protein::J => Ok('J'),
            Protein::Z => Ok('Z'),
            Protein::X => Ok('X'),
            Protein::Blank => Ok('_'),
            Protein::Pos => Ok('+'),
            Protein::Any => Ok('*'),
        }
    }

    fn str_to_vec(sequence: &str) -> Result<Vec<Protein>> {
        let mut result = Vec::<Protein>::new();

        for s in sequence.chars() {
            result.push(match Protein::match_with_char(s) {
                Ok(p) => p,
                Err(err) => return Err(err),
            });
        }

        Ok(result)
    }

    fn vec_to_str(sequence: &[Protein]) -> Result<String> {
        let mut result: Vec<char> = vec![];

        for p in sequence.iter() {
            result.push(match Protein::convert_to_char(p) {
                Ok(s) => s,
                Err(err) => return Err(err),
            });
        }

        Ok(result.into_iter().collect())
    }

    fn from_u8_vec(vec: &[u8]) -> Result<Vec<Protein>> {
        let mut result: Vec<Protein> = vec![];

        for elem in vec.iter() {
            result.push(match Protein::match_with_char(*elem as char) {
                Ok(v) => v,
                Err(err) => return Err(err),
            });
        }

        Ok(result)
    }

    fn from_u8_vec_with_freqs(vec: &[u8]) -> Result<(Vec<Protein>, Array1<f64>)> {
        let mut result: Vec<Protein> = vec![];
        let mut freqs = Array1::<f64>::zeros(Self::volume());

        for elem in vec.iter() {
            result.push(match Protein::match_with_char(*elem as char) {
                Ok(v) => {
                    freqs[v as usize] += 1f64;
                    v
                }
                Err(_) => {
                    continue;
                }
            });
        }

        freqs /= result.len() as f64;
        Ok((result, freqs))
    }

    fn from_u8_vec_with_freqs_and_indices(
        vec: &[u8],
    ) -> Result<(Vec<Protein>, Array1<f64>, Vec<Index>)> {
        let mut result: Vec<Protein> = vec![];
        let mut freqs = Array1::<f64>::zeros(Self::volume());
        let mut indices = vec![];

        let mut pass = true;
        let mut count = 0;
        let mut local_count = 0;
        for (i, elem) in vec.iter().enumerate() {
            result.push(match Protein::match_with_char(*elem as char) {
                Ok(v) => {
                    freqs[v as usize] += 1f64;
                    if !pass {
                        indices.push(Index {
                            coord: i - count,
                            offset: count,
                            local_offset: local_count,
                        });
                        local_count = 0;
                        pass = true;
                    }
                    v
                }
                Err(_) => {
                    pass = false;
                    count += 1;
                    local_count += 1;
                    continue;
                }
            });
        }

        indices.reverse();

        freqs /= result.len() as f64;
        Ok((result, freqs, indices))
    }

    fn random_seq(length: usize) -> Result<Vec<Protein>> {
        let mut thread_rng = thread_rng();

        let mut result = vec![];
        for _ in 0..length {
            result.push(thread_rng.gen_range(0..Protein::volume()).into())
        }

        Ok(result)
    }

    fn random_seq_with_freqs(length: usize) -> Result<(Vec<Protein>, Array1<f64>)> {
        let mut thread_rng = thread_rng();
        let mut freqs = Array1::<f64>::zeros(Self::volume());

        let mut result = vec![];
        for _ in 0..length {
            let elem = thread_rng.gen_range(0..Protein::volume());
            result.push(elem.into());
            freqs[elem] += 1f64;
        }

        Ok((result, freqs))
    }

    fn blank() -> Protein {
        Protein::Blank
    }

    fn pos() -> Protein {
        Protein::Pos
    }

    fn volume() -> usize {
        24
    }
}

impl BioData for DNA {
    fn match_with_char(symbol: char) -> Result<DNA> {
        match symbol {
            'A' => Ok(DNA::A),
            'T' => Ok(DNA::T),
            'C' => Ok(DNA::C),
            'G' => Ok(DNA::G),
            '_' => Ok(DNA::Blank),
            '+' => Ok(DNA::Pos),
            _ => Err(Error::CharIsNotMatchable),
        }
    }

    fn convert_to_char(dna: &DNA) -> Result<char> {
        match dna {
            DNA::A => Ok('A'),
            DNA::T => Ok('T'),
            DNA::C => Ok('C'),
            DNA::G => Ok('G'),
            DNA::Blank => Ok('_'),
            DNA::Pos => Ok('+'),
            DNA::Any => Ok('*'),
        }
    }

    fn str_to_vec(sequence: &str) -> Result<Vec<DNA>> {
        let mut result = Vec::<DNA>::new();

        for s in sequence.chars() {
            result.push(match DNA::match_with_char(s) {
                Ok(p) => p,
                Err(err) => return Err(err),
            });
        }

        Ok(result)
    }

    fn vec_to_str(sequence: &[DNA]) -> Result<String> {
        let mut result: Vec<char> = vec![];

        for p in sequence.iter() {
            result.push(match DNA::convert_to_char(p) {
                Ok(s) => s,
                Err(err) => return Err(err),
            });
        }

        Ok(result.into_iter().collect())
    }

    fn from_u8_vec(vec: &[u8]) -> Result<Vec<DNA>> {
        let mut result: Vec<DNA> = vec![];

        for elem in vec.iter() {
            result.push(match DNA::match_with_char(*elem as char) {
                Ok(v) => v,
                Err(_) => {
                    continue;
                }
            });
        }

        Ok(result)
    }

    fn from_u8_vec_with_freqs(vec: &[u8]) -> Result<(Vec<DNA>, Array1<f64>)> {
        let mut result: Vec<DNA> = vec![];
        let mut freqs = Array1::<f64>::zeros(Self::volume());

        for elem in vec.iter() {
            result.push(match DNA::match_with_char(*elem as char) {
                Ok(v) => {
                    freqs[v as usize] += 1f64;
                    v
                }
                Err(_) => {
                    continue;
                }
            });
        }

        freqs /= result.len() as f64;
        Ok((result, freqs))
    }

    fn from_u8_vec_with_freqs_and_indices(
        vec: &[u8],
    ) -> Result<(Vec<DNA>, Array1<f64>, Vec<Index>)> {
        let mut result: Vec<DNA> = vec![];
        let mut freqs = Array1::<f64>::zeros(Self::volume());
        let mut indices = vec![];

        let mut pass = true;
        let mut count = 0;
        let mut local_count = 0;
        for (i, elem) in vec.iter().enumerate() {
            result.push(match DNA::match_with_char(*elem as char) {
                Ok(v) => {
                    freqs[v as usize] += 1f64;
                    if !pass {
                        indices.push(Index {
                            coord: i - count,
                            offset: count,
                            local_offset: local_count,
                        });
                        local_count = 0;
                        pass = true;
                    }
                    v
                }
                Err(_) => {
                    pass = false;
                    count += 1;
                    local_count += 1;
                    continue;
                }
            });
        }

        indices.reverse();

        freqs /= result.len() as f64;
        Ok((result, freqs, indices))
    }

    fn random_seq(length: usize) -> Result<Vec<DNA>> {
        let mut thread_rng = thread_rng();

        let mut result = vec![];
        for _ in 0..length {
            result.push(thread_rng.gen_range(0..DNA::volume()).into())
        }

        Ok(result)
    }

    fn random_seq_with_freqs(length: usize) -> Result<(Vec<DNA>, Array1<f64>)> {
        let mut thread_rng = thread_rng();
        let mut freqs = Array1::<f64>::zeros(Self::volume());

        let mut result = vec![];
        for _ in 0..length {
            let elem = thread_rng.gen_range(0..DNA::volume());
            result.push(elem.into());
            freqs[elem] += 1f64;
        }

        Ok((result, freqs))
    }

    fn blank() -> DNA {
        DNA::Blank
    }

    fn pos() -> DNA {
        DNA::Pos
    }

    fn volume() -> usize {
        4
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Index {
    pub coord: usize,
    pub offset: usize,
    pub local_offset: usize,
}
