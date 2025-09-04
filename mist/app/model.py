import dataclasses
import json
from enum import Enum
from typing import Any

ALLELE_MISSING = '-'

class CustomEncoder(json.JSONEncoder):
    """
    Custom encoder that can handle dumping dataclass objects.
    """
    def default(self, obj: Any) -> Any:
        """
        Default encoding method.
        :param obj: Input object
        :return: Encoded object
        """
        if isinstance(obj, Enum):
            return obj.value
        if dataclasses.is_dataclass(obj):
            return dataclasses.asdict(obj)
        return super().default(obj)

class Tag(Enum):
    """
    Tags to denote the different result types.
    """
    NOVEL = 'NOVEL'
    ABSENT = 'ABSENT'
    INDEL = 'INDEL'
    EDGE = 'EDGE'


@dataclasses.dataclass(frozen=True, unsafe_hash=True)
class Profile:
    """
    Dataclass to hold ST profile data.
    """
    name: str
    metadata: list[tuple[str, str]] = dataclasses.field(hash=False)
    alleles: dict[str, str] | None = dataclasses.field(hash=False)


@dataclasses.dataclass
class Alignment:
    """
    Represents an alignment in the input sequence.
    """
    seq_id: str
    start: int
    end: int
    strand: str

    @property
    def length(self) -> int:
        """
        Returns the length of the alignment.
        :return: Alignment length
        """
        return self.end - self.start


@dataclasses.dataclass
class AlleleResult:
    """
    Match to an allele.
    """
    allele: str
    alignment: Alignment
    sequence: str | None = None
    closest_alleles: list[str] | None = None


@dataclasses.dataclass
class QueryResult:
    """
    Combined output of a locus query, can contain multiple allele matches
    """
    allele_str: str
    allele_results: list[AlleleResult]
    tags: list[Tag]
