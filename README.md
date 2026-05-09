# PubMed Author Research Tracker

A Google Sheets-based PubMed tracking pipeline that monitors research publications for a configurable list of faculty/authors, validates author identity and affiliation, removes duplicate papers, separates low-confidence matches, and generates research analytics including authorship counts, article types, collaboration indicators, and AI-generated publication themes.

This project was designed for departments, labs, and research groups that want a lightweight way to track faculty publications over time without manually searching PubMed.

---

## Features

- Searches PubMed for each tracked author using multiple name variants
- Pulls article metadata directly from PubMed XML
- Supports author aliases such as:
  - `Vishal S Kapadia`
  - `Vishal Kapadia`
  - `Kapadia VS`
  - `Kapadia V`
- Checks whether the matched author appears in the article author list or investigator list
- Checks whether article affiliations match the institution of interest
- Separates low-confidence affiliation matches into a separate sheet
- Deduplicates publications by PMID so the same paper is not counted multiple times
- Generates one master publication table
- Generates individual author tabs
- Computes analytics such as:
  - reporting period in months
  - total unique articles
  - median and range of publications per faculty member
  - percent of faculty with more than 1, 3, and 5 publications
  - first-author publication count
  - last-author publication count
  - article type counts
  - number of publications involving other institutions
  - top 5 AI-generated research themes
- Uses Google Sheets as the main user-facing interface
- Can be run manually, through Windows Task Scheduler, or through GitHub Actions

---

## Project Structure

```text
pubmed-author-tracker/
│
├── src/
│   └── script.py
│
├── logs/
│   └── debug_*.txt
│
├── out/
│   └── master.csv
│
├── requirements.txt
├── .env
└── README.md
```

---

## How It Works

The pipeline uses a Google Sheet as its configuration and output dashboard.

Users enter faculty information into a `Config` sheet. The script reads that configuration, searches PubMed, fetches article XML, filters matches, writes results back to Google Sheets, and prints analytics to stdout.

The pipeline has two date modes.

### Incremental Mode

If `starting_date` is blank, the script uses `last_run_utc` from the Google Sheet and searches from the last run date to the current date.

After a successful run, it updates `last_run_utc`.

### Manual Override Mode

If `starting_date` is filled in, the script searches from `starting_date` to `ending_date`.

In this mode, `last_run_utc` is not updated.

---

## Google Sheet Setup

Create a Google Sheet with a tab named:

```text
Config
```

The `Config` sheet should contain two sections:

```text
[SETTINGS]
key	value
starting_date	
ending_date	
last_run_utc	
ncbi_tool	pubmed-author-tracker
ncbi_email	your-email@example.com

[AUTHORS]
full_name	other_names	affiliations
Vishal S Kapadia	Vishal S Kapadia|Kapadia VS|Vishal Kapadia|Kapadia V	UT Southwestern|UTSW|University of Texas Southwestern|University of Texas Southwestern Medical Center
Myra H Wyckoff	Myra H Wyckoff|Wyckoff MH|Myra Wyckoff|Wyckoff M	UT Southwestern|UTSW|University of Texas Southwestern|University of Texas Southwestern Medical Center
```

Use the pipe character `|` to separate multiple aliases or affiliations.

---

## Output Sheets

The script writes several output tabs to the same Google Sheet.

### `Master`

Contains high-confidence articles only.

Low-confidence articles are removed from this tab.

### `Low Confidence Articles`

Contains articles where the author may match, but the affiliation is missing or does not clearly match the institution.

These are usually articles that should be manually reviewed.

### Individual Author Tabs

For each author with high-confidence publications, the script creates a separate tab using the author’s full name.

Example:

```text
Vishal S Kapadia
Myra H Wyckoff
Lina F Chalak
```

---

## Deduplication Logic

The script deduplicates articles by PMID.

This is important because the same article may appear multiple times if:

- multiple tracked faculty are authors on the same paper
- multiple aliases for the same author return the same PMID
- the author appears in both author and collaborator/investigator fields

After all article rows are collected, the dataframe is deduplicated using:

```python
df = df.drop_duplicates(subset=["pmid"])
```

This prevents the same publication from being counted multiple times in the analytics.

---

## Match Confidence Logic

The script separates matches into several status categories:

```text
AUTHOR MATCH
COLLABORATOR MATCH
AFF MATCH
AFF MISMATCH
AFF NOT FOUND
AUTHOR MISMATCH
COLLABORATOR MISMATCH
NOT FOUND
```

Articles with these statuses are treated as low confidence:

```text
AFF NOT FOUND
AFF MISMATCH
```

Those articles are sent to the `Low Confidence Articles` tab instead of the main `Master` tab.

---

## Analytics Generated

At the end of the run, the script prints a JSON analytics summary.

Example:

```json
{
  "period_months": 3,
  "total_articles": 42,
  "per_faculty": {
    "median": 1.0,
    "range": {
      "min": 0,
      "max": 6
    }
  },
  "faculty_productivity_percentages": {
    ">1_publication": 25.0,
    ">3_publications": 8.33,
    ">5_publications": 2.08
  },
  "authorship": {
    "first_author_publications": 5,
    "last_author_publications": 12
  },
  "collaboration": {
    "articles_with_other_institutions": 18
  },
  "article_types": {
    "Journal Article": 30,
    "Review": 5,
    "Clinical Trial": 2
  },
  "themes": {
    "themes": [
      {
        "rank": 1,
        "name": "Neonatal outcomes",
        "article_count": 12,
        "description": "Studies focused on health outcomes in newborn and neonatal populations."
      }
    ]
  }
}
```

---

## Environment Variables

Create a `.env` file in the project root.

```env
SPREADSHEET_ID=your_google_sheet_id_here
GOOGLE_APPLICATION_CREDENTIALS_JSON={"type":"service_account",...}
NCBI_API_KEY=your_ncbi_api_key_here
GEMINI_API_KEY=your_gemini_api_key_here
```

### Required

```text
SPREADSHEET_ID
GOOGLE_APPLICATION_CREDENTIALS_JSON
```

### Recommended

```text
NCBI_API_KEY
GEMINI_API_KEY
```

`NCBI_API_KEY` is recommended because it increases PubMed request limits.

`GEMINI_API_KEY` is used for AI-generated research theme extraction.

---

## Google Service Account Setup

1. Create a Google Cloud project.
2. Enable:
   - Google Sheets API
   - Google Drive API
3. Create a service account.
4. Generate a JSON key for the service account.
5. Copy the full JSON into the `GOOGLE_APPLICATION_CREDENTIALS_JSON` environment variable.
6. Share your Google Sheet with the service account email.

The service account email usually looks like:

```text
something@project-name.iam.gserviceaccount.com
```

Give it editor access to the Google Sheet.

---

## Installation

Install dependencies:

```bash
pip install -r requirements.txt
```

Example `requirements.txt`:

```txt
requests
pandas
PyYAML
tenacity
gspread
google-auth
python-dotenv
google-genai
pydantic
```

---

## Running the Script

From the project root:

```bash
python src/script.py
```

The script will:

1. Read the Google Sheet config.
2. Search PubMed.
3. Fetch article XML.
4. Match authors and affiliations.
5. Deduplicate publications.
6. Write results to Google Sheets.
7. Print analytics to stdout.
8. Update `last_run_utc` if running in incremental mode.

---

## Scheduling Quarterly Runs

### Option 1: Windows Task Scheduler

This is the simplest local option.

Create a scheduled task that runs:

```bash
python C:\path\to\pubmed-author-tracker\src\script.py
```

Set the schedule to run every 3 months.

Make sure the task starts in the project folder so `.env` loads correctly.

### Option 2: GitHub Actions

This is better for cloud automation.

A scheduled GitHub Actions workflow can run the script quarterly using repository secrets for the environment variables.

---

## Google Apps Script Email Updates

A Google Apps Script can be attached to the output Google Sheet to email a summary after the analytics cell is updated.

This allows the Python script to focus on data collection while Google Sheets handles email notifications.

---

## Current Limitations

- Author matching is rule-based and depends on aliases provided in the Config sheet.
- Common names may still require manual review.
- Affiliation matching is based on keyword matching.
- Low-confidence articles should be reviewed manually.
- AI theme extraction depends on the Gemini API and may occasionally require retry logic.
- The same PMID is counted once globally, so multi-faculty coauthored papers are not double-counted in total article counts.

---

## Future Improvements

Potential future features include:

- A Streamlit or web dashboard
- More advanced author disambiguation
- ORCID support
- Better affiliation-to-author backtracking
- Automated email reports
- Historical trend charts
- Quarterly PDF reports
- Institution-level collaboration network graphs
- Manual review buttons for low-confidence articles
- Config template auto-generation

---

## Intended Use

This project is intended for research groups, departments, or faculty teams that want an easy way to monitor PubMed publication activity over time using Google Sheets as the main dashboard.

It is especially useful when tracking multiple faculty members across recurring reporting periods such as monthly, quarterly, or yearly publication reviews.
