# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Emit ``llms.txt`` and ``llms-full.txt`` beside the built HTML.

Two files, for two ways of arriving:

``llms.txt``
    The index described at https://llmstxt.org -- a title, a summary, and a
    flat list of every page with a one-line description. Small enough to read
    whole before deciding what to fetch.

``llms-full.txt``
    The documentation as one plain-text file, **including the output of every
    executed code cell**. That last part is the reason this is generated from
    the doctree rather than shipped as a copy of the Markdown sources: the
    sources show code with nothing coming out of it, and half of what these
    pages demonstrate is what the code prints -- the fitted numbers, the error
    messages, the array that came back in the wrong order. A reader who only
    sees the input learns the API but not the behaviour.

Generated at build time for the same reason the guide pages execute: a
hand-written summary of the documentation is a second source of truth, and it
goes stale silently.
"""

from __future__ import annotations

import os

from docutils import nodes
from sphinx.util import logging

logger = logging.getLogger(__name__)


def _redact(text):
    """
    Take the build machine's home directory out of the text.

    A stray warning from a badly installed dependency prints an absolute path,
    and cell output is captured verbatim. On CI that path is meaningless; built
    locally it is somebody's home directory, in a file whose whole purpose is
    to be fetched and read elsewhere.
    """
    home = os.path.expanduser("~")
    return text.replace(home, "~") if home not in ("", "/") else text

#: Nodes whose text is already carried by a parent, or which are furniture.
_SKIP = (nodes.system_message, nodes.comment, nodes.problematic)


class _TextRenderer(nodes.NodeVisitor):
    """Doctree to plain text, keeping headings, code and code output."""

    def __init__(self, document):
        super().__init__(document)
        self.parts: list[str] = []
        self.depth = 0
        self._suppress = 0

    # -- helpers ----------------------------------------------------------

    def _emit(self, text):
        if text and not self._suppress:
            self.parts.append(text)

    def unknown_visit(self, node):
        if isinstance(node, _SKIP):
            raise nodes.SkipNode

    def unknown_departure(self, node):
        pass

    # -- structure --------------------------------------------------------

    def visit_section(self, node):
        self.depth += 1

    def depart_section(self, node):
        self.depth -= 1

    def visit_title(self, node):
        self._emit(f"\n{'#' * max(self.depth, 1)} {node.astext()}\n")
        raise nodes.SkipNode

    def visit_paragraph(self, node):
        self._emit(node.astext() + "\n")
        raise nodes.SkipNode

    def visit_literal_block(self, node):
        language = node.get('language', '')
        language = '' if language in ('default', 'none') else language
        self._emit(f"```{language}\n{node.astext()}\n```\n")
        raise nodes.SkipNode

    def visit_bullet_list(self, node):
        for item in node.children:
            self._emit(f"- {' '.join(item.astext().split())}")
        self._emit("")
        raise nodes.SkipNode

    def visit_enumerated_list(self, node):
        for number, item in enumerate(node.children, 1):
            self._emit(f"{number}. {' '.join(item.astext().split())}")
        self._emit("")
        raise nodes.SkipNode

    def visit_table(self, node):
        self._emit(node.astext() + "\n")
        raise nodes.SkipNode

    def visit_image(self, node):
        # A figure is the one thing this format cannot carry. Say so, rather
        # than dropping it and leaving a paragraph referring to nothing.
        self._emit(f"[figure: {node.get('alt') or node.get('uri', 'plot')}]\n")
        raise nodes.SkipNode

    # -- executed cell output ---------------------------------------------

    def visit_container(self, node):
        classes = node.get('classes', [])
        if 'cell_output' in classes:
            text = node.astext().strip()
            if text:
                self._emit(f"Output:\n```\n{text}\n```\n")
            raise nodes.SkipNode


def _render(app, docname):
    """One page as plain text, or ``None`` if it cannot be rendered."""
    try:
        doctree = app.env.get_and_resolve_doctree(docname, app.builder)
    except Exception as error:                               # noqa: BLE001
        logger.warning("llms.txt: could not resolve %s: %s", docname, error)
        return None
    renderer = _TextRenderer(doctree)
    doctree.walkabout(renderer)
    return "\n".join(renderer.parts).strip()


def _summary(app, docname):
    """The page's first real paragraph, for the index line."""
    try:
        doctree = app.env.get_doctree(docname)
    except Exception:                                        # noqa: BLE001
        return ""
    for paragraph in doctree.findall(nodes.paragraph):
        text = " ".join(paragraph.astext().split())
        if len(text) > 30:
            return text if len(text) <= 200 else text[:197].rsplit(' ', 1)[0] + "..."
    return ""


def _title(app, docname):
    node = app.env.titles.get(docname)
    return node.astext() if node is not None else docname


def _ordered_docnames(app):
    """Every page, in the order the toctree presents them."""
    seen, order = set(), []

    def walk(docname):
        if docname in seen:
            return
        seen.add(docname)
        order.append(docname)
        for child in app.env.toctree_includes.get(docname, []):
            walk(child)

    walk(app.config.root_doc)
    order.extend(sorted(set(app.env.all_docs) - seen))
    return order


def write_llms_txt(app, exception):
    if exception is not None or app.builder.name != 'html':
        return

    base = (app.config.html_baseurl or "").rstrip('/')
    docnames = _ordered_docnames(app)

    index = [
        f"# {app.config.project}",
        "",
        f"> {app.config.llms_txt_summary}",
        "",
        app.config.llms_txt_notes.strip(),
        "",
        "## Documentation",
        "",
    ]
    for docname in docnames:
        url = f"{base}/{docname}.html" if base else f"{docname}.html"
        summary = _summary(app, docname)
        index.append(f"- [{_title(app, docname)}]({url})"
                     + (f": {summary}" if summary else ""))
    index += [
        "",
        "## Optional",
        "",
        f"- [Everything above as one file]({base}/llms-full.txt): the same "
        f"pages in full, including the printed output of every executed "
        f"example.",
        "",
    ]

    full = [
        f"# {app.config.project} {app.config.release} -- full documentation",
        "",
        f"> {app.config.llms_txt_summary}",
        "",
        "Generated from the built documentation, so the code examples below "
        "are shown with the output they actually produced.",
        "",
    ]
    for docname in docnames:
        rendered = _render(app, docname)
        if rendered:
            full += [f"\n\n{'=' * 70}\nPage: {docname}\n{'=' * 70}\n",
                     rendered]

    out = app.outdir
    (out / "llms.txt").write_text(_redact("\n".join(index)), encoding='utf-8')
    (out / "llms-full.txt").write_text(_redact("\n".join(full)) + "\n",
                                       encoding='utf-8')
    logger.info("wrote llms.txt (%d pages) and llms-full.txt (%d KB)",
                len(docnames),
                len("\n".join(full).encode('utf-8')) // 1024)


def setup(app):
    app.add_config_value('llms_txt_summary', '', 'html')
    app.add_config_value('llms_txt_notes', '', 'html')
    app.connect('build-finished', write_llms_txt)
    return {'version': '1.0', 'parallel_read_safe': True,
            'parallel_write_safe': True}
