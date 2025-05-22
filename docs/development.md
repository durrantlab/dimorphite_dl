# Development

This comprehensive guide provides detailed instructions to help maintainers effectively develop, test, document, build, and release new versions of `dimorphite_dl`.

## Setting up the Development Environment

`dimorphite_dl` utilizes [`pixi`](https://pixi.sh/latest/) for managing environments and dependencies, streamlining the setup process. Follow these precise steps to configure your development environment:

1.  **Clone the repository:**
    Begin by obtaining a local copy of the `dimorphite_dl` codebase:

    ```bash
    git clone git@GitHub.com:durrantlab/dimorphite_dl.git
    cd dimorphite_dl
    ```
2.  **Install dependencies:**
    Install all necessary dependencies by running:

    ```bash
    pixi install
    ```
3.  **Activate the development environment:**
    To enter the isolated virtual environment configured specifically for `dimorphite_dl` development, execute:

    ```bash
    pixi shell
    ```

You are now fully prepared and equipped to develop `dimorphite_dl`.

## Code Formatting and Style Guide

Maintaining consistent style and formatting across the codebase is crucial for readability and maintainability.
`dimorphite_dl` employs automated formatting tools configured to enforce standardized style guidelines.
Execute the following command to apply formatting automatically:

```bash
pixi run format
```

This command sequentially runs `black` for Python formatting, `isort` for managing imports, and `markdownlint-cli2` to enforce markdown formatting standards, ensuring your contributions align with project conventions.

## Documentation

`dimorphite_dl`'s documentation is built using MkDocs, allowing easy creation and maintenance of high-quality documentation.
To locally preview documentation changes, serve the documentation by running:

```bash
pixi run -e docs serve-docs
```

After execution, open your web browser and visit [`http://127.0.0.1:8000/`](http://127.0.0.1:8000/) to review changes in real-time.

## Testing

Writing and maintaining tests is essential for ensuring code correctness, reliability, and stability.
Execute `dimorphite_dl`'s tests with:

```bash
pixi run -e dev tests
```

Additionally, you can evaluate test coverage to identify untested areas and improve overall reliability by running:

```bash
pixi run -e dev coverage
```

Review the generated coverage reports to address any gaps in testing.

## Building the Package

Prepare `dimorphite_dl` for publishing or distribution by building the package.
Execute:

```bash
pixi run build
```

Upon completion, inspect the `dist` directory for the generated distribution files, which are ready for publication.

## Bumping Version

Releasing a new version of `dimorphite_dl` requires updating version information, documenting changes, and creating a corresponding release tag. Follow these steps precisely to ensure consistency and traceability:

1. **Update version identifiers:**
   Change the version number in all of the following files:

   - `pyproject.toml`
   - `pixi.toml`
   - `dimorphite_dl/__init__.py` (update the `__version__` string)

   Ensure that the version is consistent across all files.

2. **Update the changelog:**
   Document all notable changes since the previous release in the `CHANGELOG.md` file. Follow a consistent and clear format to help users understand what has changed.

3. **Commit the changes:**
   Stage and commit the version bump and changelog update using a clear and standardized message, for example:

   ```bash
   git add .
   git commit -m "bump: v1.2.5"
   ```

4. **Tag the commit:**
   Create a version tag that follows the `v<version>` format:

   ```bash
   git tag v1.2.5
   git push origin main --tags
   ```

5. **Create a GitHub release:**
   Navigate to the [GitHub Releases](https://github.com/durrantlab/dimorphite_dl/releases) page and draft a new release:

   - Tag version: `v1.2.5`
   - Release title: `v1.2.5`
   - Description: Copy the relevant changelog section or summarize the key changes.

   Attach the built package files from the `dist/` directory, if desired.

## Publishing to PyPI

Once the version number is updated and the package is built, it can be published to PyPI.
Execute:

```bash
pixi run publish
```

For preliminary testing or release candidates, it is highly recommended to publish to TestPyPI first.
Execute:

```bash
pixi run publish-test
```

Publishing to TestPyPI allows you to validate packaging correctness and installation processes without affecting production users.

## Maintenance Best Practices

To maintain high quality and reliability of `dimorphite_dl`, adhere to the following best practices:

-   Regularly synchronize your local repository with the main branch to incorporate the latest updates:

    ```bash
    git pull origin main
    ```
-   Frequently review and address open issues and pull requests on GitHub.
-   Clearly document changes in commit messages, issue descriptions, and pull requests.
-   Routinely verify dependencies and update them as necessary to maintain compatibility and security.

Adhering to these guidelines ensures a robust, stable, and continuously improving `dimorphite_dl` project.

This expanded documentation guide covers the entire workflow comprehensively, providing clarity and precision for effective `dimorphite_dl` project maintenance.
