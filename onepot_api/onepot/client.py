import json
from collections.abc import Generator
from typing import Literal, NotRequired, TypedDict

import httpx


class BBFilter(TypedDict):
    """Per-(reaction_class, bb_index) Tanimoto similarity bound on candidate BBs.

    Use in `bb_filters` on `Client.search()` / `Client.search_stream()`. The
    `reaction_class` value is the `rxn_<hash>` string from a prior
    `decompose=True` response — pass it through verbatim. Either
    similarity bound may be omitted; omit both to leave the position
    unconstrained.
    """

    reaction_class: str
    bb_index: int
    min_similarity: NotRequired[float]
    max_similarity: NotRequired[float]


class PropertyRange(TypedDict, total=False):
    """Inclusive min/max bounds for one screening-index property."""

    min: float
    max: float


class SpacePropertyFilters(TypedDict, total=False):
    """Property filters accepted by :meth:`Client.sample_space`."""

    molecular_weight: PropertyRange
    clogp: PropertyRange
    tpsa: PropertyRange
    hbd: PropertyRange
    hba: PropertyRange
    rotatable_bonds: PropertyRange
    heavy_atoms: PropertyRange
    fraction_csp3: PropertyRange
    aromatic_rings: PropertyRange
    rings: PropertyRange
    qed: PropertyRange


class SpaceMoleculeProperties(TypedDict):
    molecular_weight: float
    clogp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    heavy_atoms: int
    fraction_csp3: float
    aromatic_rings: int
    rings: int
    qed: float
    murcko_scaffold: str
    generic_murcko_scaffold: str


class SpaceMolecule(TypedDict):
    smiles: str
    inchikey: str
    price_usd: Literal[125, 295]
    properties: NotRequired[SpaceMoleculeProperties]


class SpaceSampleResponse(TypedDict):
    molecules: list[SpaceMolecule]
    seed: int


class Client:
    def __init__(self, api_key: str, base_url: str = "https://api.onepot.ai"):
        self.api_key = api_key
        self.base_url = base_url.rstrip("/")
        self.client = httpx.Client(
            headers={"X-API-Key": api_key},
            timeout=240.0,
            follow_redirects=True,
        )

    def search(
        self,
        smiles_list: list[str],
        max_results: int = 100,
        substructure_search: bool = False,
        max_depth: int = 1,
        include_chemistry_risk: bool = False,
        include_chemistry_risk_score: bool = False,
        max_price: int | None = None,
        max_supplier_risk: str | None = None,
        max_chemistry_risk: str | None = None,
        decompose: bool = False,
        bb_filters: list[BBFilter] | None = None,
        exact_lookup: bool = False,
    ):
        """Search for purchasable analogs of the given molecules.

        Args:
            smiles_list: List of SMILES strings to search.
            max_results: Maximum results per query (default 100).
            substructure_search: Match substructures instead of similarity.
            max_depth: Synthesis depth (currently only 1 supported).
            include_chemistry_risk: Include chemistry_risk field (low/medium/high).
            include_chemistry_risk_score: Include raw probability score (0-1).
            max_price: Maximum price in USD. Results above this are excluded.
            max_supplier_risk: Maximum supplier risk level ("low", "medium", or "high").
            max_chemistry_risk: Maximum chemistry risk level ("low", "medium", or "high").
            decompose: If True, the response includes a `decompositions` block per
                query listing every retrosynthetic path considered, with the
                `reaction_class` and the BB SMILES at each position.
                Pricing unchanged.
            bb_filters: Per-(reaction_class, bb_index) Tanimoto bounds applied to
                candidate BBs before enumeration. Each entry is a `BBFilter`
                TypedDict: {reaction_class, bb_index, min_similarity?,
                max_similarity?}. Pricing unchanged.
            exact_lookup: If True, price each query molecule directly — a catalog
                hit or the cheapest single-step decomposition — and skip the
                analog/similarity search entirely. Each query returns at most one
                result: the query molecule itself (similarity 1.0), or none if it
                can't be priced. Much faster for bulk pricing of pre-enumerated
                libraries. Cannot be combined with substructure_search,
                decompose, or bb_filters.

        Returns:
            Dict with:
                - queries: List of query results, each containing query_smiles,
                  query_inchikey, and results (list of {smiles, inchikey, similarity,
                  supplier_risk, price_usd, and optionally chemistry_risk /
                  chemistry_risk_score}). Enumerated results are also tagged with
                  `reaction_class` and `bbs` (list of {bb_index, smiles}) when the
                  search synthesized them. When `decompose=True`, each query entry
                  also has a `decompositions` field.
                - credits_used: Number of credits consumed.
                - credits_remaining: Remaining credit balance.
        """
        body = {
            "smiles_list": smiles_list,
            "max_results": max_results,
            "substructure_search": substructure_search,
            "max_depth": max_depth,
            "include_chemistry_risk": include_chemistry_risk,
            "include_chemistry_risk_score": include_chemistry_risk_score,
            "decompose": decompose,
            "exact_lookup": exact_lookup,
        }
        if max_price is not None:
            body["max_price"] = max_price
        if max_supplier_risk is not None:
            body["max_supplier_risk"] = max_supplier_risk
        if max_chemistry_risk is not None:
            body["max_chemistry_risk"] = max_chemistry_risk
        if bb_filters is not None:
            body["bb_filters"] = bb_filters
        resp = self.client.post(f"{self.base_url}/v1/search", json=body)
        if not resp.is_success:
            print(f"Error {resp.status_code}: {resp.text}")
        resp.raise_for_status()
        return resp.json()

    def search_stream(
        self,
        smiles: str,
        max_results: int = 100,
        substructure_search: bool = False,
        include_chemistry_risk: bool = False,
        include_chemistry_risk_score: bool = False,
        max_price: int | None = None,
        max_supplier_risk: str | None = None,
        max_chemistry_risk: str | None = None,
        decompose: bool = False,
        bb_filters: list[BBFilter] | None = None,
        exact_lookup: bool = False,
    ) -> Generator[dict]:
        """Stream search results for a single molecule via Server-Sent Events.

        Yields dicts with 'status' and 'message' keys as the search progresses.
        Status lifecycle: starting → synthesis → rescoring → complete.
        The final event (status='complete') includes a 'results' list with the
        same fields as the batch search endpoint. When `decompose=True` it also
        includes a `decompositions` list.

        Args:
            smiles: Single SMILES string to search.
            max_results: Maximum results (default 100).
            substructure_search: Match substructures instead of similarity.
            include_chemistry_risk: Include chemistry_risk field (low/medium/high).
            include_chemistry_risk_score: Include raw probability score (0-1).
            max_price: Maximum price in USD. Results above this are excluded.
            max_supplier_risk: Maximum supplier risk level ("low", "medium", or "high").
            max_chemistry_risk: Maximum chemistry risk level ("low", "medium", or "high").
            decompose: If True, the final event includes a `decompositions` list
                describing the retrosynthetic paths considered, each with a
                `reaction_class` and the BB SMILES at each position.
                Pricing unchanged.
            bb_filters: Per-(reaction_class, bb_index) Tanimoto bounds applied to
                candidate BBs before enumeration. Each entry is a `BBFilter`
                TypedDict: {reaction_class, bb_index, min_similarity?,
                max_similarity?}. Pricing unchanged.
            exact_lookup: If True, price the query molecule directly — a catalog
                hit or the cheapest single-step decomposition — and skip the
                analog/similarity search. The final event's 'results' holds at
                most one entry: the query molecule itself (similarity 1.0), or
                none if it can't be priced. The status lifecycle is
                starting → pricing → complete. Cannot be combined with
                substructure_search, decompose, or bb_filters.

        Yields:
            Dict with 'status' and 'message'. Final event also contains 'results'
            (and `decompositions` when decompose=True). Enumerated results are
            tagged with `reaction_class` and `bbs` (list of {bb_index, smiles}).
        """
        body = {
            "smiles": smiles,
            "max_results": max_results,
            "substructure_search": substructure_search,
            "include_chemistry_risk": include_chemistry_risk,
            "include_chemistry_risk_score": include_chemistry_risk_score,
            "decompose": decompose,
            "exact_lookup": exact_lookup,
        }
        if max_price is not None:
            body["max_price"] = max_price
        if max_supplier_risk is not None:
            body["max_supplier_risk"] = max_supplier_risk
        if max_chemistry_risk is not None:
            body["max_chemistry_risk"] = max_chemistry_risk
        if bb_filters is not None:
            body["bb_filters"] = bb_filters
        with self.client.stream(
            "POST",
            f"{self.base_url}/v1/search/stream",
            json=body,
        ) as resp:
            if not resp.is_success:
                raise httpx.HTTPStatusError(f"Error {resp.status_code}", request=resp.request, response=resp)
            buffer = ""
            for chunk in resp.iter_text():
                buffer += chunk
                while "\n\n" in buffer:
                    event, buffer = buffer.split("\n\n", 1)
                    for line in event.strip().split("\n"):
                        if line.startswith("data: "):
                            yield json.loads(line[6:])

    def order(
        self,
        smiles: list[str],
        email: str,
        notes: str = "",
    ):
        """Submit an order for the given molecules.

        Args:
            smiles: List of SMILES strings to order.
            email: Contact email for the order.
            notes: Optional free-text notes.

        Returns:
            Dict with:
                - order_id: UUID string identifying the order.
                - molecule_count: Number of molecules in the order.
        """
        resp = self.client.post(
            f"{self.base_url}/v1/order",
            json={
                "smiles": smiles,
                "email": email,
                "notes": notes,
            },
        )
        if not resp.is_success:
            print(f"Error {resp.status_code}: {resp.text}")
        resp.raise_for_status()
        return resp.json()

    def sample_space(
        self,
        count: int = 384,
        strategy: Literal["diverse", "random"] = "random",
        seed: int | None = None,
        properties: SpacePropertyFilters | None = None,
        include_properties: bool = False,
        exclude_inchikeys: list[str] | None = None,
        exclude_generic_scaffolds: list[str] | None = None,
    ) -> SpaceSampleResponse:
        """Sample makeable molecules for a screening campaign.

        Args:
            count: Molecules to return, from 1 to 10,000 (default 384).
            strategy: ``"diverse"`` for generic-Murcko-balanced sampling or
                ``"random"`` for a seeded random sample of the indexed pool.
            seed: Optional uint32 seed. Omit to have the API generate and return
                one; reuse a returned seed to reproduce the sample.
            properties: Optional inclusive ranges for molecular_weight, clogp,
                tpsa, hbd, hba, rotatable_bonds, heavy_atoms, fraction_csp3,
                aromatic_rings, rings, and qed.
            include_properties: Include descriptors and Murcko frameworks on
                each returned molecule (default False).
            exclude_inchikeys: Optional identifiers to omit from a follow-up
                batch (maximum 100,000).
            exclude_generic_scaffolds: Optional exact generic Bemis-Murcko
                scaffold SMILES to omit (maximum 10,000). Values can be copied
                from prior molecule ``properties`` blocks.

        Returns:
            Dict containing exactly ``count`` molecules and the effective
            ``seed``. When no seed is supplied, store the returned value to
            replay the request against the same serving index.

        Raises:
            httpx.HTTPStatusError: If authentication fails, request validation
                fails, or the server cannot complete the request. Validation
                errors include the server's detail message.
        """
        body = {
            "count": count,
            "strategy": strategy,
            "include_properties": include_properties,
        }
        if seed is not None:
            body["seed"] = seed
        if properties is not None:
            body["properties"] = properties
        if exclude_inchikeys is not None:
            body["exclude_inchikeys"] = exclude_inchikeys
        if exclude_generic_scaffolds is not None:
            body["exclude_generic_scaffolds"] = exclude_generic_scaffolds

        resp = self.client.post(f"{self.base_url}/v1/space/sample", json=body)
        if not resp.is_success:
            try:
                detail = resp.json().get("detail", resp.text)
            except (ValueError, AttributeError):
                detail = resp.text
            raise httpx.HTTPStatusError(
                f"Error {resp.status_code}: {detail}",
                request=resp.request,
                response=resp,
            )
        return resp.json()
