"""Data loaders for the de novo sequencing task with optional ion mobility."""

import functools
import os
from typing import List, Optional, Tuple

import lightning.pytorch as pl
import numpy as np
import torch
from depthcharge.data import AnnotatedSpectrumIndex

from ..data.datasets import (
    AnnotatedSpectrumDataset,
    SpectrumDataset,
    IMAnnotatedSpectrumDataset,
    IMSpectrumDataset,
)


class DeNovoDataModule(pl.LightningDataModule):
    """
    Data loader to prepare MS/MS spectra for a Spec2Pep predictor.

    Parameters
    ----------
    train_index : Optional[AnnotatedSpectrumIndex]
        The spectrum index file corresponding to the training data.
    valid_index : Optional[AnnotatedSpectrumIndex]
        The spectrum index file corresponding to the validation data.
    test_index : Optional[AnnotatedSpectrumIndex]
        The spectrum index file corresponding to the testing data.
    train_batch_size : int
        The batch size to use for training.
    eval_batch_size : int
        The batch size to use for inference.
    n_peaks : Optional[int]
        The number of top-n most intense peaks to keep in each spectrum. `None`
        retains all peaks.
    min_mz : float
        The minimum m/z to include. The default is 140 m/z, in order to exclude
        TMT and iTRAQ reporter ions.
    max_mz : float
        The maximum m/z to include.
    min_intensity : float
        Remove peaks whose intensity is below `min_intensity` percentage of the
        base peak intensity.
    remove_precursor_tol : float
        Remove peaks within the given mass tolerance in Dalton around the
        precursor mass.
    n_workers : int, optional
        The number of workers to use for data loading. By default, the number of
        available CPU cores on the current machine is used.
    random_state : Optional[int]
        The NumPy random state. ``None`` leaves mass spectra in the order they
        were parsed.
    use_ion_mobility : bool
        Whether to load and pass ion mobility values through the data pipeline.
        Default False for backward compatibility.
    im_array_train : Optional[np.ndarray]
        Ion mobility values for training spectra.
    im_array_valid : Optional[np.ndarray]
        Ion mobility values for validation spectra.
    im_array_test : Optional[np.ndarray]
        Ion mobility values for test spectra.
    """

    def __init__(
        self,
        train_index: Optional[AnnotatedSpectrumIndex] = None,
        valid_index: Optional[AnnotatedSpectrumIndex] = None,
        test_index: Optional[AnnotatedSpectrumIndex] = None,
        train_batch_size: int = 128,
        eval_batch_size: int = 1028,
        n_peaks: Optional[int] = 150,
        min_mz: float = 50.0,
        max_mz: float = 2500.0,
        min_intensity: float = 0.01,
        remove_precursor_tol: float = 2.0,
        n_workers: Optional[int] = None,
        random_state: Optional[int] = None,
        use_ion_mobility: bool = False,
        im_array_train: Optional[np.ndarray] = None,
        im_array_valid: Optional[np.ndarray] = None,
        im_array_test: Optional[np.ndarray] = None,
    ):
        super().__init__()
        self.train_index = train_index
        self.valid_index = valid_index
        self.test_index = test_index
        self.train_batch_size = train_batch_size
        self.eval_batch_size = eval_batch_size
        self.n_peaks = n_peaks
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.min_intensity = min_intensity
        self.remove_precursor_tol = remove_precursor_tol
        self.n_workers = n_workers if n_workers is not None else os.cpu_count()
        self.rng = np.random.default_rng(random_state)
        self.train_dataset = None
        self.valid_dataset = None
        self.test_dataset = None
        self.use_ion_mobility = use_ion_mobility
        self.im_array_train = im_array_train
        self.im_array_valid = im_array_valid
        self.im_array_test = im_array_test

    def setup(self, stage: str = None, annotated: bool = True) -> None:
        """
        Set up the PyTorch Datasets.

        Parameters
        ----------
        stage : str {"fit", "validate", "test"}
            The stage indicating which Datasets to prepare. All are prepared by
            default.
        annotated: bool
            True if peptide sequence annotations are available for the test
            data.
        """
        if stage in (None, "fit", "validate"):
            if self.use_ion_mobility:
                make_dataset = functools.partial(
                    IMAnnotatedSpectrumDataset,
                    n_peaks=self.n_peaks,
                    min_mz=self.min_mz,
                    max_mz=self.max_mz,
                    min_intensity=self.min_intensity,
                    remove_precursor_tol=self.remove_precursor_tol,
                )
                if self.train_index is not None:
                    self.train_dataset = make_dataset(
                        self.train_index,
                        im_array=self.im_array_train,
                        random_state=self.rng,
                    )
                if self.valid_index is not None:
                    self.valid_dataset = make_dataset(
                        self.valid_index,
                        im_array=self.im_array_valid,
                    )
            else:
                make_dataset = functools.partial(
                    AnnotatedSpectrumDataset,
                    n_peaks=self.n_peaks,
                    min_mz=self.min_mz,
                    max_mz=self.max_mz,
                    min_intensity=self.min_intensity,
                    remove_precursor_tol=self.remove_precursor_tol,
                )
                if self.train_index is not None:
                    self.train_dataset = make_dataset(
                        self.train_index,
                        random_state=self.rng,
                    )
                if self.valid_index is not None:
                    self.valid_dataset = make_dataset(self.valid_index)

        if stage in (None, "test"):
            if self.use_ion_mobility:
                if annotated:
                    make_dataset = functools.partial(
                        IMAnnotatedSpectrumDataset,
                        n_peaks=self.n_peaks,
                        min_mz=self.min_mz,
                        max_mz=self.max_mz,
                        min_intensity=self.min_intensity,
                        remove_precursor_tol=self.remove_precursor_tol,
                    )
                else:
                    make_dataset = functools.partial(
                        IMSpectrumDataset,
                        n_peaks=self.n_peaks,
                        min_mz=self.min_mz,
                        max_mz=self.max_mz,
                        min_intensity=self.min_intensity,
                        remove_precursor_tol=self.remove_precursor_tol,
                    )
                if self.test_index is not None:
                    self.test_dataset = make_dataset(
                        self.test_index,
                        im_array=self.im_array_test,
                    )
            else:
                make_dataset = functools.partial(
                    AnnotatedSpectrumDataset if annotated else SpectrumDataset,
                    n_peaks=self.n_peaks,
                    min_mz=self.min_mz,
                    max_mz=self.max_mz,
                    min_intensity=self.min_intensity,
                    remove_precursor_tol=self.remove_precursor_tol,
                )
                if self.test_index is not None:
                    self.test_dataset = make_dataset(self.test_index)

    def _make_loader(
        self,
        dataset: torch.utils.data.Dataset,
        batch_size: int,
        shuffle: bool = False,
    ) -> torch.utils.data.DataLoader:
        """
        Create a PyTorch DataLoader.

        Parameters
        ----------
        dataset : torch.utils.data.Dataset
            A PyTorch Dataset.
        batch_size : int
            The batch size to use.
        shuffle : bool
            Option to shuffle the batches.

        Returns
        -------
        torch.utils.data.DataLoader
            A PyTorch DataLoader.
        """
        collate_fn = (
            prepare_batch_im if self.use_ion_mobility else prepare_batch
        )
        return torch.utils.data.DataLoader(
            dataset,
            batch_size=batch_size,
            collate_fn=collate_fn,
            pin_memory=True,
            num_workers=self.n_workers,
            shuffle=shuffle,
        )

    def train_dataloader(self) -> torch.utils.data.DataLoader:
        """Get the training DataLoader."""
        return self._make_loader(
            self.train_dataset, self.train_batch_size, shuffle=True
        )

    def val_dataloader(self) -> torch.utils.data.DataLoader:
        """Get the validation DataLoader."""
        return self._make_loader(self.valid_dataset, self.eval_batch_size)

    def test_dataloader(self) -> torch.utils.data.DataLoader:
        """Get the test DataLoader."""
        return self._make_loader(self.test_dataset, self.eval_batch_size)

    def predict_dataloader(self) -> torch.utils.data.DataLoader:
        """Get the predict DataLoader."""
        return self._make_loader(self.test_dataset, self.eval_batch_size)


def prepare_batch(
    batch: List[Tuple[torch.Tensor, float, int, str]]
) -> Tuple[torch.Tensor, torch.Tensor, np.ndarray]:
    """
    Collate MS/MS spectra into a batch (original 3-column precursors).

    Parameters
    ----------
    batch : List[Tuple[torch.Tensor, float, int, str]]
        A batch of data from an AnnotatedSpectrumDataset, consisting of for each
        spectrum (i) a tensor with the m/z and intensity peak values, (ii), the
        precursor m/z, (iii) the precursor charge, (iv) the spectrum identifier.

    Returns
    -------
    spectra : torch.Tensor of shape (batch_size, n_peaks, 2)
        The padded mass spectra tensor with the m/z and intensity peak values
        for each spectrum.
    precursors : torch.Tensor of shape (batch_size, 3)
        A tensor with the precursor neutral mass, precursor charge, and
        precursor m/z.
    spectrum_ids : np.ndarray
        The spectrum identifiers (during de novo sequencing) or peptide
        sequences (during training).
    """
    spectra, precursor_mzs, precursor_charges, spectrum_ids = list(zip(*batch))
    spectra = torch.nn.utils.rnn.pad_sequence(spectra, batch_first=True)
    precursor_mzs = torch.tensor(precursor_mzs)
    precursor_charges = torch.tensor(precursor_charges)
    precursor_masses = (precursor_mzs - 1.007276) * precursor_charges
    precursors = torch.vstack(
        [precursor_masses, precursor_charges, precursor_mzs]
    ).T.float()
    return spectra, precursors, np.asarray(spectrum_ids)


def prepare_batch_im(
    batch: List[Tuple[torch.Tensor, float, int, float, str]]
) -> Tuple[torch.Tensor, torch.Tensor, np.ndarray]:
    """
    Collate MS/MS spectra into a batch with ion mobility (4-column precursors).

    The batch tuples have 5 elements: (spectrum, precursor_mz, precursor_charge,
    ion_mobility, spectrum_id/peptide).

    Parameters
    ----------
    batch : List[Tuple[torch.Tensor, float, int, float, str]]
        A batch of data from an IMSpectrumDataset or IMAnnotatedSpectrumDataset.

    Returns
    -------
    spectra : torch.Tensor of shape (batch_size, n_peaks, 2)
        The padded mass spectra tensor.
    precursors : torch.Tensor of shape (batch_size, 4)
        A tensor with (precursor_mass, charge, m/z, 1/K0).
    spectrum_ids : np.ndarray
        The spectrum identifiers or peptide sequences.
    """
    spectra, precursor_mzs, precursor_charges, ion_mobilities, spectrum_ids = (
        list(zip(*batch))
    )
    spectra = torch.nn.utils.rnn.pad_sequence(spectra, batch_first=True)
    precursor_mzs = torch.tensor(precursor_mzs)
    precursor_charges = torch.tensor(precursor_charges)
    ion_mobilities = torch.tensor(ion_mobilities)
    precursor_masses = (precursor_mzs - 1.007276) * precursor_charges
    precursors = torch.vstack(
        [precursor_masses, precursor_charges, precursor_mzs, ion_mobilities]
    ).T.float()
    return spectra, precursors, np.asarray(spectrum_ids)
