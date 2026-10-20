"""
Body site normalization utilities for pathogen classification.

Provides standardized body site names based on the pathogen sheet categories
and normalizes various synonyms and related terms.
"""

def normalize_body_site(site):
    """
    Normalize body site names to standard categories used in pathogen classification.

    Args:
        site (str): Raw body site name (e.g., "gut", "stool", "CSF", "plasma")

    Returns:
        str: Normalized body site name in lowercase

    Standard Categories (from pathogen_sheet.csv pathogenic_sites):
        - abscess
        - blood
        - bone
        - brain
        - cornea
        - csf
        - ear
        - gut
        - lung
        - lymph
        - oral
        - resp (respiratory)
        - skin
        - sputum
        - stool
        - urine
        - urogenital
        - sterile (special case - any normally sterile site)
    """

    if not site:
        return "unknown"

    # Convert to lowercase and strip whitespace
    site = str(site).lower().strip()

    # Define normalization mappings
    SITE_MAPPINGS = {
        # Blood and blood products
        'blood': ['blood', 'serum', 'plasma', 'erythrocyte', 'leukocyte',
                  'whole blood', 'peripheral blood', 'bloodstream', 'bacteremia',
                  'sepsis', 'septicemia', 'hematogenous'],

        # Gastrointestinal - Gut
        'gut': ['gut', 'intestine', 'intestinal', 'gi', 'gastrointestinal',
                'bowel', 'colon', 'colonic', 'duodenum', 'jejunum', 'ileum',
                'cecum', 'digestive', 'enteric', 'gastric', 'stomach'],

        # Gastrointestinal - Stool/Feces
        'stool': ['stool', 'feces', 'fecal', 'faeces', 'faecal', 'diarrhea',
                  'diarrhoea', 'rectal', 'rectum'],

        # Urinary System
        'urine': ['urine', 'urinary', 'bladder', 'urinary bladder', 'uti',
                  'urinary tract', 'nephric', 'renal pelvis', 'cystitis'],

        # Urogenital (broader than just urine)
        'urogenital': ['urogenital', 'genitourinary', 'gu', 'genital',
                       'reproductive', 'vaginal', 'vagina', 'cervical', 'cervix',
                       'uterine', 'uterus', 'prostatic', 'prostate', 'testicular',
                       'epididymal', 'seminal'],

        # Respiratory - General
        'resp': ['resp', 'respiratory', 'respir', 'pulmonary', 'bronchial',
                 'bronchi', 'tracheal', 'trachea', 'airway', 'airways'],

        # Respiratory - Lung specific
        'lung': ['lung', 'lungs', 'pulmonary', 'pneumonia', 'pneumonic',
                 'alveolar', 'pleural', 'pleura', 'pleurisy'],

        # Respiratory - Sputum
        'sputum': ['sputum', 'phlegm', 'mucus', 'expectoration', 'bronchial secretion'],

        # Central Nervous System
        'csf': ['csf', 'cerebrospinal fluid', 'cerebrospinal', 'spinal fluid',
                'meningeal', 'meninges', 'meningitis'],

        'brain': ['brain', 'cerebral', 'cerebrum', 'encephalic', 'encephalitis',
                  'intracranial', 'neural', 'neurological'],

        # Skin and Soft Tissue
        'skin': ['skin', 'cutaneous', 'dermal', 'dermis', 'epidermal', 'epidermis',
                 'subcutaneous', 'integument', 'wound', 'cellulitis'],

        # Abscesses and Pus
        'abscess': ['abscess', 'pus', 'purulent', 'pustule', 'furuncle',
                    'carbuncle', 'empyema', 'pyogenic'],

        # Oral/Mouth
        'oral': ['oral', 'mouth', 'buccal', 'gingival', 'gum', 'gums',
                 'dental', 'tooth', 'teeth', 'periodontal', 'pharyngeal',
                 'pharynx', 'throat', 'tonsillar', 'tonsil', 'tonsils',
                 'oropharyngeal', 'nasopharyngeal'],

        # Ear
        'ear': ['ear', 'otic', 'aural', 'otitis', 'mastoid'],

        # Eye
        'cornea': ['cornea', 'corneal', 'eye', 'ocular', 'ophthalmic',
                   'conjunctival', 'conjunctiva', 'keratitis'],

        # Lymphatic System
        'lymph': ['lymph', 'lymphatic', 'lymph node', 'lymph nodes',
                  'lymphadenitis', 'adenitis'],

        # Bone and Joints
        'bone': ['bone', 'osseous', 'skeletal', 'osteomyelitis',
                 'joint', 'articular', 'arthritis', 'synovial', 'synovium'],

        # Nasal
        'nasal': ['nasal', 'nose', 'nares', 'nostril', 'sinus', 'sinuses',
                  'paranasal', 'rhinitis', 'nasopharynx'],

        # Sterile Sites (normally sterile body sites)
        'sterile': ['sterile', 'normally sterile', 'sterile site', 'sterile fluid',
                    'clean', 'peritoneal', 'peritoneum', 'ascites', 'pleural fluid',
                    'pericardial', 'pericardium', 'synovial fluid'],
    }

    # Check each mapping
    for standard_site, synonyms in SITE_MAPPINGS.items():
        # First check exact match
        if site in synonyms:
            return standard_site

    # Then check partial matches (but be more careful)
    for standard_site, synonyms in SITE_MAPPINGS.items():
        for synonym in synonyms:
            # Only match if the site contains the full synonym as a word
            # This prevents "meningitis" matching "gut" because it contains "giti"
            if len(synonym) >= 4 and synonym in site:
                return standard_site

    # If no match found, return the original (normalized to lowercase)
    return site


def get_pathogen_classification(pathogen_info, sample_site):
    """
    Determine pathogen classification based on sample site context.

    Args:
        pathogen_info (dict): Pathogen information with fields:
            - callclass: general classification
            - pathogenic_sites: list of sites where pathogenic
            - commensal_sites: list of sites where commensal
        sample_site (str): The sample collection site (will be normalized)

    Returns:
        tuple: (classification, is_direct_match)
            - classification: "primary", "opportunistic", "potential", "commensal", or "unknown"
            - is_direct_match: True if sample_site matches pathogenic/commensal sites

    Logic:
        1. If sample_site == "sterile" → use callclass (all pathogens pathogenic in sterile sites)
        2. If sample_site in pathogenic_sites → pathogenic (use callclass)
        3. If sample_site in commensal_sites → "commensal"
        4. Otherwise → use callclass but mark as not direct match
    """

    # Normalize the sample site
    normalized_site = normalize_body_site(sample_site)
    # Get pathogen data
    callclass = (pathogen_info.get('callclass') or pathogen_info.get('general_classification') or '').lower().strip()
    pathogenic_sites = pathogen_info.get('pathogenic_sites', [])
    commensal_sites = pathogen_info.get('commensal_sites', [])

    # Normalize pathogenic and commensal sites
    if isinstance(pathogenic_sites, str):
        pathogenic_sites = [s.strip() for s in pathogenic_sites.split(',')]
    if isinstance(commensal_sites, str):
        commensal_sites = [s.strip() for s in commensal_sites.split(',')]

    pathogenic_sites = [normalize_body_site(s) for s in pathogenic_sites if s]
    commensal_sites = [normalize_body_site(s) for s in commensal_sites if s]

    # Special case: sterile sites
    if normalized_site == 'sterile':
        # Everything is pathogenic in sterile sites
        return callclass or 'primary', True
    # Check if sample site is in pathogenic sites
    # Also check if the normalized site matches any normalized pathogenic site
    if normalized_site in pathogenic_sites:
        # Use the general classification
        return 'primary', True

    # Check if sample site is in commensal sites
    # Also check if the normalized site matches any normalized commensal site
    if normalized_site in commensal_sites:
        return 'commensal', True

    # Special handling for gut/stool equivalence
    # If sample is stool and gut is in commensal sites, it's commensal
    if normalized_site == 'stool' and 'gut' in commensal_sites:
        return 'commensal', True
    if normalized_site == 'gut' and 'stool' in commensal_sites:
        return 'commensal', True

    # Not a direct match - use general classification but mark as indirect
    return callclass or 'unknown', False


# Test cases
if __name__ == "__main__":
    print("Testing body site normalization:")
    print("-" * 60)

    test_sites = [
        "blood", "plasma", "serum", "bacteremia",
        "gut", "stool", "feces", "intestine", "GI",
        "urine", "UTI", "bladder",
        "urogenital", "vaginal", "prostate",
        "resp", "lung", "pulmonary", "sputum",
        "CSF", "meningitis", "brain",
        "skin", "wound", "cellulitis",
        "oral", "throat", "dental",
        "nasal", "nose", "sinus",
        "sterile", "peritoneal",
        "Unknown Site", "bone", "abscess"
    ]

    for site in test_sites:
        normalized = normalize_body_site(site)
        print(f"{site:20s} → {normalized}")

    print("\n" + "=" * 60)
    print("Testing pathogen classification:")
    print("-" * 60)

    # Example pathogen
    test_pathogen = {
        'callclass': 'opportunistic',
        'pathogenic_sites': 'blood, resp, urine',
        'commensal_sites': 'gut, skin'
    }

    test_samples = [
        ('blood', 'Should be opportunistic (pathogenic site)'),
        ('plasma', 'Should be opportunistic (blood synonym)'),
        ('gut', 'Should be commensal (commensal site)'),
        ('stool', 'Should be commensal (gut synonym)'),
        ('oral', 'Should be opportunistic (not listed, use callclass)'),
        ('sterile', 'Should be opportunistic (sterile site rule)'),
    ]

    for sample_site, description in test_samples:
        classification, direct = get_pathogen_classification(test_pathogen, sample_site)
        match_type = "DIRECT" if direct else "INDIRECT"
        print(f"{sample_site:15s} → {classification:15s} [{match_type}] - {description}")
