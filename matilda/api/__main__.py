"""Start the API server: ``python -m matilda.api --host 127.0.0.1 --port 8000``.

Limits are read from MATILDA_API_* environment variables, see
:class:`matilda.api.app.Settings`.
"""

import argparse


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="python -m matilda.api",
        description="Serve matilda's harmonic persistent homology over HTTP.",
    )
    parser.add_argument("--host", default="127.0.0.1", help="by default 127.0.0.1")
    parser.add_argument("--port", type=int, default=8000, help="by default 8000")
    args = parser.parse_args(argv)

    import uvicorn

    from .app import create_app

    uvicorn.run(create_app(), host=args.host, port=args.port)


if __name__ == "__main__":
    main()
