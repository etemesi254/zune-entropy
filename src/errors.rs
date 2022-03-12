use std::error::Error;
use std::fmt::{Debug, Display, Formatter};

#[derive(Clone)]
pub enum EntropyErrors
{
    IOError(String),
    CorruptHeader(String),
    CorruptStream(String),
}

impl Debug for EntropyErrors
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result
    {
        match self
        {
            Self::IOError(ref r) => writeln!(f, "Input/Output Error,Reason:{}", r),
            Self::CorruptStream(ref r) => writeln!(f, "Corrupt stream Reason: {}", r),
            Self::CorruptHeader(ref r) => writeln!(f, "Corrupt headers, Reason :{}", r),
        }
    }
}

impl Display for EntropyErrors
{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result
    {
        match self
        {
            Self::IOError(ref r) => writeln!(f, "Input/Output Error,Reason:{}", r),
            Self::CorruptStream(ref r) => writeln!(f, "Corrupt stream Reason: {}", r),
            Self::CorruptHeader(ref r) => writeln!(f, "Corrupt headers, Reason :{}", r),
        }
    }
}

impl Error for EntropyErrors {}

impl From<std::io::Error> for EntropyErrors
{
    fn from(err: std::io::Error) -> Self
    {
        Self::IOError(err.to_string())
    }
}
