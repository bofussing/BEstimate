import typing as t
import importlib.resources
import gzip
from pathlib import Path
import io

import pandas as pd

from BEstimate import constants

if t.TYPE_CHECKING:
    from importlib.resources.abc import Traversable

_MODULE_NAME = __name__


class DataFiles:
    MD5_H_SAPIENS_INTERFACES_GZIP: t.ClassVar[str] = (
        constants.CHECKSUM_MD5_H_SAPIENS_INTERFACES_GZIP
    )
    # Checksum from http://interactomeinsider.yulab.org/downloads/interfacesHQ/H_sapiens_interfacesHQ.txt
    # as of 2025-08-17
    MD5_H_SAPIENS_INTERFACES_TXT: t.ClassVar[str] = (
        constants.CHECKSUM_MD5_H_SAPIENS_INTERFACES_TXT
    )
    H_SAPIENS_INTERFACES_COLUMNS: t.ClassVar[tuple[str, ...]] = (
        "P1",
        "P2",
        "Source",
        "P1_IRES",
        "P2_IRES",
    )

    @classmethod
    def get_home_sapiens_interfaces_path(cls) -> Path:
        """Get the path to the h_sapiens_interfaces.txt.gz file.

        The providence of this file is from
        http://interactomeinsider.yulab.org/downloads.html for 'Highest
        Confidence Interfaces'. Specifically, the file
        http://interactomeinsider.yulab.org/downloads/interfacesHQ/H_sapiens_interfacesHQ.txt

        :return: h_sapiens_interfaces.txt.gz file path
        """
        as_resource: "Traversable" = importlib.resources.files(_MODULE_NAME).joinpath(
            "H_sapiens_interfaces.txt.gz"
        )
        as_path = Path(str(as_resource))
        return as_path

    @classmethod
    def get_home_sapiens_interfaces_as_dataframe(cls) -> "pd.DataFrame":
        """Decompress and parse the h_sapiens_interfaces.txt.gz file as a pandas DataFrame.

        The providence of this file is from
        http://interactomeinsider.yulab.org/downloads.html for 'Highest
        Confidence Interfaces'. Specifically, the file
        http://interactomeinsider.yulab.org/downloads/interfacesHQ/H_sapiens_interfacesHQ.txt

        :return: DataFrame containing the h_sapiens_interfaces data
        """
        file_path = cls.get_home_sapiens_interfaces_path()
        with gzip.open(file_path, mode="rb") as f:
            file_as_bytes = f.read()
            with io.StringIO(file_as_bytes.decode("utf-8")) as text_f:
                df = pd.read_table(text_f, sep="\t")
        return df
