#!/usr/bin/env python3
"""Create a Nature-style broadness figure from AlphaGenome strict-skin box stats."""

from __future__ import annotations

import argparse
import csv
import html
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


WIDTH = 2100
HEIGHT = 1450
SQUARE_SIZE = 1500
FONT_FAMILY = "Arial, Helvetica, DejaVu Sans, sans-serif"

BLACK = "#222222"
GRAY = "#8D8D8D"
PALE_GRAY = "#F1F1F1"
TEAL = "#009E73"
ORANGE = "#D55E00"
SOFT_LAVENDER = "#B9B4C7"
DEEP_BLUE = "#168C8C"


@dataclass(frozen=True)
class SnvStats:
    input_snv: str
    n: float
    minimum: float
    q1: float
    median: float
    q3: float
    maximum: float

    @property
    def iqr(self) -> float:
        return self.q3 - self.q1

    @property
    def full_range(self) -> float:
        return self.maximum - self.minimum


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render a Nature-style figure showing how broad AlphaGenome RNA "
            "strict-skin scores are within individual SNVs."
        )
    )
    parser.add_argument(
        "--stats-csv",
        required=True,
        type=Path,
        help="alphagenome_strict_skin_boxplot_stats.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory. Defaults to <stats-csv parent>/nature_visualizations.",
    )
    parser.add_argument(
        "--basename",
        default="alphagenome_strict_skin_broadness",
        help="Output filename stem.",
    )
    return parser.parse_args()


def parse_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def load_stats(csv_path: Path) -> tuple[list[dict[str, str]], list[SnvStats]]:
    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        raw_rows = list(reader)

    required = [
        "input_snv",
        "rna_abs_quantile_n",
        "rna_abs_quantile_min",
        "rna_abs_quantile_q1",
        "rna_abs_quantile_median",
        "rna_abs_quantile_q3",
        "rna_abs_quantile_max",
    ]
    missing = [column for column in required if column not in (reader.fieldnames or [])]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    records: list[SnvStats] = []
    for row in raw_rows:
        values = {
            "n": parse_float(row.get("rna_abs_quantile_n")),
            "minimum": parse_float(row.get("rna_abs_quantile_min")),
            "q1": parse_float(row.get("rna_abs_quantile_q1")),
            "median": parse_float(row.get("rna_abs_quantile_median")),
            "q3": parse_float(row.get("rna_abs_quantile_q3")),
            "maximum": parse_float(row.get("rna_abs_quantile_max")),
        }
        if any(value is None for value in values.values()):
            continue
        records.append(
            SnvStats(
                input_snv=row["input_snv"],
                n=float(values["n"]),
                minimum=float(values["minimum"]),
                q1=float(values["q1"]),
                median=float(values["median"]),
                q3=float(values["q3"]),
                maximum=float(values["maximum"]),
            )
        )
    return raw_rows, records


def quantile(values: Iterable[float], q: float) -> float:
    sorted_values = sorted(values)
    if not sorted_values:
        return float("nan")
    position = q * (len(sorted_values) - 1)
    low = int(math.floor(position))
    high = int(math.ceil(position))
    if low == high:
        return sorted_values[low]
    weight = position - low
    return sorted_values[low] * (1 - weight) + sorted_values[high] * weight


def fmt_float(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def short_snv_label(value: str) -> str:
    if len(value) <= 16:
        return value
    chrom, _, rest = value.partition(":")
    return f"{chrom}:{rest[:6]}..."


def pct(fraction: float) -> str:
    return f"{100 * fraction:.1f}%"


def hex_to_rgb(color: str, alpha: float = 1.0) -> tuple[int, int, int, int]:
    color = color.lstrip("#")
    return (
        int(color[0:2], 16),
        int(color[2:4], 16),
        int(color[4:6], 16),
        int(255 * alpha),
    )


class SvgCanvas:
    def __init__(self, width: int, height: int) -> None:
        self.width = width
        self.height = height
        width_mm = 183
        height_mm = width_mm * height / width
        self.parts: list[str] = [
            (
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_mm:.1f}mm" height="{height_mm:.1f}mm" '
                f'viewBox="0 0 {width} {height}">'
            ),
            '<rect width="100%" height="100%" fill="white"/>',
        ]

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        color: str = BLACK,
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        self.parts.append(
            (
                f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
                f'stroke="{color}" stroke-width="{width:.2f}" stroke-opacity="{alpha:.3f}" '
                'stroke-linecap="round"/>'
            )
        )

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        fill: str = "none",
        stroke: str = "none",
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        self.parts.append(
            (
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" '
                f'fill="{fill}" fill-opacity="{alpha:.3f}" stroke="{stroke}" '
                f'stroke-width="{width:.2f}"/>'
            )
        )

    def circle(self, x: float, y: float, r: float, fill: str, alpha: float = 1.0) -> None:
        self.parts.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{fill}" fill-opacity="{alpha:.3f}"/>'
        )

    def polyline(
        self,
        points: list[tuple[float, float]],
        color: str,
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        point_text = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        self.parts.append(
            (
                f'<polyline points="{point_text}" fill="none" stroke="{color}" '
                f'stroke-width="{width:.2f}" stroke-opacity="{alpha:.3f}" stroke-linejoin="round"/>'
            )
        )

    def path(
        self,
        commands: list[str],
        color: str,
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        self.parts.append(
            (
                f'<path d="{" ".join(commands)}" fill="none" stroke="{color}" '
                f'stroke-width="{width:.2f}" stroke-opacity="{alpha:.3f}" stroke-linecap="round"/>'
            )
        )

    def text(
        self,
        x: float,
        y: float,
        text: str,
        size: float,
        color: str = BLACK,
        weight: str = "normal",
        anchor: str = "start",
    ) -> None:
        self.parts.append(
            (
                f'<text x="{x:.2f}" y="{y:.2f}" fill="{color}" font-family="{FONT_FAMILY}" '
                f'font-size="{size:.1f}" font-weight="{weight}" text-anchor="{anchor}">'
                f"{html.escape(text)}</text>"
            )
        )

    def save(self, path: Path) -> None:
        self.parts.append("</svg>")
        path.write_text("\n".join(self.parts), encoding="utf-8")


class PngCanvas:
    def __init__(self, width: int, height: int) -> None:
        from PIL import Image, ImageDraw, ImageFont

        self.Image = Image
        self.ImageDraw = ImageDraw
        self.ImageFont = ImageFont
        self.width = width
        self.height = height
        self.image = Image.new("RGBA", (width, height), "white")
        self.draw = ImageDraw.Draw(self.image)
        self.font_cache: dict[tuple[float, str], object] = {}

    def font(self, size: float, weight: str = "normal"):
        key = (size, weight)
        if key in self.font_cache:
            return self.font_cache[key]
        candidates = []
        if weight == "bold":
            candidates.extend(
                [
                    "/path/to/fonts/arialbd.ttf",
                    "/path/to/fonts/calibrib.ttf",
                    "/path/to/fonts/DejaVuSans-Bold.ttf",
                ]
            )
        candidates.extend(
            [
                "/path/to/fonts/arial.ttf",
                "/path/to/fonts/calibri.ttf",
                "/path/to/fonts/DejaVuSans.ttf",
            ]
        )
        font_obj = None
        for candidate in candidates:
            if Path(candidate).exists():
                font_obj = self.ImageFont.truetype(candidate, int(round(size)))
                break
        if font_obj is None:
            font_obj = self.ImageFont.load_default()
        self.font_cache[key] = font_obj
        return font_obj

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        color: str = BLACK,
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        self.draw.line(
            [(x1, y1), (x2, y2)],
            fill=hex_to_rgb(color, alpha),
            width=max(1, int(round(width))),
        )

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        fill: str = "none",
        stroke: str = "none",
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        fill_color = None if fill == "none" else hex_to_rgb(fill, alpha)
        stroke_color = None if stroke == "none" else hex_to_rgb(stroke, 1.0)
        if fill_color is not None and alpha < 1.0:
            overlay = self.Image.new("RGBA", (self.width, self.height), (255, 255, 255, 0))
            overlay_draw = self.ImageDraw.Draw(overlay)
            overlay_draw.rectangle([x, y, x + w, y + h], fill=fill_color)
            self.image = self.Image.alpha_composite(self.image, overlay)
            self.draw = self.ImageDraw.Draw(self.image)
            fill_color = None
        self.draw.rectangle(
            [x, y, x + w, y + h],
            fill=fill_color,
            outline=stroke_color,
            width=max(1, int(round(width))),
        )

    def circle(self, x: float, y: float, r: float, fill: str, alpha: float = 1.0) -> None:
        self.draw.ellipse([x - r, y - r, x + r, y + r], fill=hex_to_rgb(fill, alpha))

    def polyline(
        self,
        points: list[tuple[float, float]],
        color: str,
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        self.draw.line(points, fill=hex_to_rgb(color, alpha), width=max(1, int(round(width))))

    def path(
        self,
        commands: list[str],
        color: str,
        width: float = 1.0,
        alpha: float = 1.0,
    ) -> None:
        segments: list[tuple[float, float, float, float]] = []
        current: tuple[float, float] | None = None
        for command in commands:
            tokens = command.split()
            if len(tokens) != 3:
                continue
            _, x_text, y_text = tokens
            point = (float(x_text), float(y_text))
            if tokens[0] == "M":
                current = point
            elif tokens[0] == "L" and current is not None:
                segments.append((current[0], current[1], point[0], point[1]))
        for x1, y1, x2, y2 in segments:
            self.line(x1, y1, x2, y2, color=color, width=width, alpha=alpha)

    def text(
        self,
        x: float,
        y: float,
        text: str,
        size: float,
        color: str = BLACK,
        weight: str = "normal",
        anchor: str = "start",
    ) -> None:
        font = self.font(size, weight)
        bbox = self.draw.textbbox((0, 0), text, font=font)
        text_width = bbox[2] - bbox[0]
        if anchor == "middle":
            x -= text_width / 2
        elif anchor == "end":
            x -= text_width
        self.draw.text((x, y - size), text, fill=hex_to_rgb(color, 1.0), font=font)

    def save(self, path: Path) -> None:
        self.image.convert("RGB").save(path, dpi=(300, 300))


def x_scale(x0: float, width: float, xmin: float = 0.0, xmax: float = 1.0):
    return lambda value: x0 + (value - xmin) / (xmax - xmin) * width


def y_scale(y0: float, height: float, ymin: float = 0.0, ymax: float = 1.0):
    return lambda value: y0 + height - (value - ymin) / (ymax - ymin) * height


def draw_axes(
    canvas,
    x0: float,
    y0: float,
    width: float,
    height: float,
    xticks: list[float],
    yticks: list[float],
    xlabel: str,
    ylabel: str,
    xformatter=lambda x: f"{x:g}",
    yformatter=lambda y: f"{y:g}",
    tick_size: float = 20,
    label_size: float = 24,
) -> tuple:
    sx = x_scale(x0, width, min(xticks), max(xticks))
    sy = y_scale(y0, height, min(yticks), max(yticks))
    for tick in yticks:
        y = sy(tick)
        canvas.line(x0, y, x0 + width, y, color=PALE_GRAY, width=1.0)
        canvas.line(x0 - 8, y, x0, y, color=BLACK, width=1.2)
        canvas.text(x0 - 16, y + 7, yformatter(tick), tick_size, color=BLACK, anchor="end")
    for tick in xticks:
        x = sx(tick)
        canvas.line(x, y0 + height, x, y0 + height + 8, color=BLACK, width=1.2)
        canvas.text(x, y0 + height + 34, xformatter(tick), tick_size, color=BLACK, anchor="middle")
    canvas.line(x0, y0, x0, y0 + height, color=BLACK, width=1.5)
    canvas.line(x0, y0 + height, x0 + width, y0 + height, color=BLACK, width=1.5)
    canvas.text(x0 + width / 2, y0 + height + 72, xlabel, label_size, color=BLACK, anchor="middle")
    canvas.text(x0, y0 - 18, ylabel, label_size, color=BLACK, anchor="start")
    return sx, sy


def draw_panel_label(canvas, x: float, y: float, label: str, title: str) -> None:
    canvas.text(x, y, label, 34, weight="bold")
    canvas.text(x + 48, y, title, 28, weight="bold")


def draw_ranked_intervals(canvas, records: list[SnvStats]) -> None:
    panel_x, panel_y, panel_w, panel_h = 150, 205, 1250, 500
    draw_panel_label(canvas, panel_x - 35, panel_y - 68, "a", "RNA scores per SNV span most of the quantile scale")
    xticks = [0, 1000, 2000, 3000, 4000, len(records)]
    sx, sy = draw_axes(
        canvas,
        panel_x,
        panel_y,
        panel_w,
        panel_h,
        xticks=xticks,
        yticks=[0, 0.25, 0.5, 0.75, 1.0],
        xlabel="SNVs in input CSV order",
        ylabel="RNA absolute quantile",
        xformatter=lambda x: f"{int(x):,}",
        yformatter=lambda y: f"{y:.2f}",
    )
    minmax_commands = []
    iqr_commands = []
    for idx, record in enumerate(records):
        x = sx(idx)
        minmax_commands.extend([f"M {x:.2f} {sy(record.minimum):.2f}", f"L {x:.2f} {sy(record.maximum):.2f}"])
        iqr_commands.extend([f"M {x:.2f} {sy(record.q1):.2f}", f"L {x:.2f} {sy(record.q3):.2f}"])
    canvas.path(minmax_commands, color=GRAY, width=1.0, alpha=0.22)
    canvas.path(iqr_commands, color=TEAL, width=1.8, alpha=0.62)
    canvas.line(panel_x + panel_w - 260, panel_y + 38, panel_x + panel_w - 205, panel_y + 38, color=GRAY, width=4, alpha=0.35)
    canvas.text(panel_x + panel_w - 190, panel_y + 47, "min-max", 20, color=BLACK)
    canvas.line(panel_x + panel_w - 260, panel_y + 74, panel_x + panel_w - 205, panel_y + 74, color=TEAL, width=4, alpha=0.8)
    canvas.text(panel_x + panel_w - 190, panel_y + 83, "IQR", 20, color=BLACK)


def draw_ranked_intervals_square(canvas, records: list[SnvStats]) -> None:
    panel_x, panel_y, panel_w, panel_h = 170, 210, 1180, 920
    canvas.text(92, 105, "a", 38, weight="bold")
    canvas.text(145, 105, "RNA scores per SNV span most of the quantile scale", 34, weight="bold")
    sx, sy = draw_axes(
        canvas,
        panel_x,
        panel_y,
        panel_w,
        panel_h,
        xticks=[0, 1000, 2000, 3000, 4000, len(records)],
        yticks=[0, 0.25, 0.5, 0.75, 1.0],
        xlabel="SNVs in input CSV order",
        ylabel="",
        xformatter=lambda x: f"{int(x):,}",
        yformatter=lambda y: f"{y:.2f}",
        tick_size=26,
        label_size=32,
    )
    minmax_commands = []
    iqr_commands = []
    for idx, record in enumerate(records):
        x = sx(idx)
        minmax_commands.extend([f"M {x:.2f} {sy(record.minimum):.2f}", f"L {x:.2f} {sy(record.maximum):.2f}"])
        iqr_commands.extend([f"M {x:.2f} {sy(record.q1):.2f}", f"L {x:.2f} {sy(record.q3):.2f}"])

    canvas.path(minmax_commands, color=SOFT_LAVENDER, width=1.1, alpha=0.40)
    canvas.path(iqr_commands, color=DEEP_BLUE, width=1.7, alpha=0.54)

    legend_x = panel_x + panel_w - 280
    legend_y = panel_y + 58
    canvas.rect(
        legend_x - 34,
        legend_y - 38,
        284,
        108,
        fill="#FFFFFF",
        stroke="#D8D2DE",
        width=1.2,
        alpha=0.78,
    )
    canvas.line(legend_x, legend_y, legend_x + 62, legend_y, color=SOFT_LAVENDER, width=5, alpha=0.75)
    canvas.text(legend_x + 82, legend_y + 8, "min-max", 22, color=BLACK)
    canvas.line(legend_x, legend_y + 42, legend_x + 62, legend_y + 42, color=DEEP_BLUE, width=5, alpha=0.85)
    canvas.text(legend_x + 82, legend_y + 50, "IQR", 22, color=BLACK)


def draw_ecdf(canvas, records: list[SnvStats]) -> None:
    panel_x, panel_y, panel_w, panel_h = 1550, 205, 420, 500
    draw_panel_label(canvas, panel_x - 35, panel_y - 68, "b", "Within-SNV width")
    sx, sy = draw_axes(
        canvas,
        panel_x,
        panel_y,
        panel_w,
        panel_h,
        xticks=[0, 0.25, 0.5, 0.75, 1.0],
        yticks=[0, 0.25, 0.5, 0.75, 1.0],
        xlabel="Width on quantile scale",
        ylabel="Fraction of SNVs",
        xformatter=lambda x: f"{x:.2f}",
        yformatter=lambda y: f"{y:.2f}",
    )
    for values, color, label, y_offset in [
        ([record.iqr for record in records], TEAL, "IQR", 0),
        ([record.full_range for record in records], ORANGE, "min-max", 36),
    ]:
        sorted_values = sorted(values)
        points = [(sx(value), sy((idx + 1) / len(sorted_values))) for idx, value in enumerate(sorted_values)]
        canvas.polyline(points, color=color, width=3.2, alpha=0.92)
        canvas.line(panel_x + 35, panel_y + 42 + y_offset, panel_x + 88, panel_y + 42 + y_offset, color=color, width=4)
        canvas.text(panel_x + 102, panel_y + 51 + y_offset, label, 20, color=BLACK)


def draw_representative_boxes(canvas, records: list[SnvStats]) -> None:
    panel_x, panel_y, panel_w, panel_h = 150, 900, 1180, 390
    draw_panel_label(canvas, panel_x - 35, panel_y - 68, "c", "Representative per-SNV box summaries")
    axis_x = panel_x + 315
    sx, sy = draw_axes(
        canvas,
        axis_x,
        panel_y,
        panel_w - 315,
        panel_h,
        xticks=[0, 0.25, 0.5, 0.75, 1.0],
        yticks=[0, 1],
        xlabel="RNA absolute quantile",
        ylabel="",
        xformatter=lambda x: f"{x:.2f}",
        yformatter=lambda y: "",
    )
    sorted_records = sorted(records, key=lambda item: item.iqr)
    quantile_labels = [(0.05, "P05"), (0.25, "P25"), (0.50, "P50"), (0.75, "P75"), (0.95, "P95"), (0.99, "P99")]
    picks: list[tuple[str, SnvStats]] = []
    for q, label in quantile_labels:
        idx = int(round(q * (len(sorted_records) - 1)))
        picks.append((label, sorted_records[idx]))

    top = panel_y + 36
    row_gap = 54
    for idx, (label, record) in enumerate(picks):
        y = top + idx * row_gap
        canvas.text(panel_x, y + 7, f"{label} IQR", 18, color=BLACK)
        canvas.text(panel_x + 105, y + 7, short_snv_label(record.input_snv), 18, color=BLACK)
        canvas.line(sx(record.minimum), y, sx(record.maximum), y, color=GRAY, width=2.2, alpha=0.7)
        canvas.rect(sx(record.q1), y - 15, sx(record.q3) - sx(record.q1), 30, fill=TEAL, stroke=TEAL, alpha=0.34)
        canvas.line(sx(record.median), y - 18, sx(record.median), y + 18, color=BLACK, width=2.2)
        canvas.line(sx(record.minimum), y - 8, sx(record.minimum), y + 8, color=GRAY, width=2.0)
        canvas.line(sx(record.maximum), y - 8, sx(record.maximum), y + 8, color=GRAY, width=2.0)
        canvas.text(sx(1.0) + 15, y + 7, f"IQR={record.iqr:.3f}", 18, color=GRAY)


def draw_takeaway_panel(canvas, raw_rows: list[dict[str, str]], records: list[SnvStats]) -> None:
    panel_x, panel_y, panel_w, panel_h = 1550, 900, 420, 390
    draw_panel_label(canvas, panel_x - 35, panel_y - 68, "d", "Key structure and spread")

    unique_snvs = len({row.get("input_snv", "") for row in raw_rows if row.get("input_snv", "")})
    rna_n = [record.n for record in records]
    iqr_values = [record.iqr for record in records]
    range_values = [record.full_range for record in records]
    iqr_gt_030 = sum(value > 0.30 for value in iqr_values) / len(iqr_values)
    range_gt_090 = sum(value > 0.90 for value in range_values) / len(range_values)

    atac_n_values = {
        parse_float(row.get("atac_abs_quantile_n"))
        for row in raw_rows
        if parse_float(row.get("atac_abs_quantile_n")) is not None
    }
    atac_n_label = "1 per SNV" if atac_n_values == {1.0} else f"{len(atac_n_values)} values"

    items = [
        ("Unique SNVs", f"{unique_snvs:,}"),
        ("RNA SNVs with summaries", f"{len(records):,}"),
        ("RNA outputs per SNV", f"median {quantile(rna_n, 0.50):.0f}"),
        ("Median RNA IQR", fmt_float(quantile(iqr_values, 0.50))),
        ("SNVs with IQR > 0.30", pct(iqr_gt_030)),
        ("Median min-max range", fmt_float(quantile(range_values, 0.50))),
        ("SNVs with range > 0.90", pct(range_gt_090)),
        ("ATAC after strict skin", atac_n_label),
    ]
    y = panel_y + 15
    for label, value in items:
        canvas.text(panel_x, y + 24, label, 22, color=BLACK)
        canvas.text(panel_x + panel_w, y + 24, value, 24, color=BLACK, weight="bold", anchor="end")
        canvas.line(panel_x, y + 40, panel_x + panel_w, y + 40, color=PALE_GRAY, width=1.2)
        y += 45


def draw_figure(canvas, raw_rows: list[dict[str, str]], records: list[SnvStats]) -> None:
    canvas.text(
        WIDTH / 2,
        52,
        "AlphaGenome strict-skin RNA results remain broad within individual SNVs",
        32,
        weight="bold",
        anchor="middle",
    )
    canvas.text(
        WIDTH / 2,
        88,
        "Each interval summarizes the absolute quantile scores returned for one SNV across retained RNA genes/tracks",
        22,
        color=GRAY,
        anchor="middle",
    )
    draw_ranked_intervals(canvas, records)
    draw_ecdf(canvas, records)
    draw_representative_boxes(canvas, records)
    draw_takeaway_panel(canvas, raw_rows, records)


def write_svg(path: Path, raw_rows: list[dict[str, str]], records: list[SnvStats]) -> None:
    canvas = SvgCanvas(WIDTH, HEIGHT)
    draw_figure(canvas, raw_rows, records)
    canvas.save(path)


def write_png(path: Path, raw_rows: list[dict[str, str]], records: list[SnvStats]) -> bool:
    try:
        canvas = PngCanvas(WIDTH, HEIGHT)
    except ModuleNotFoundError:
        return False
    draw_figure(canvas, raw_rows, records)
    canvas.save(path)
    return True


def write_panel_a_square_svg(path: Path, records: list[SnvStats]) -> None:
    canvas = SvgCanvas(SQUARE_SIZE, SQUARE_SIZE)
    draw_ranked_intervals_square(canvas, records)
    canvas.save(path)


def write_panel_a_square_png(path: Path, records: list[SnvStats]) -> bool:
    try:
        canvas = PngCanvas(SQUARE_SIZE, SQUARE_SIZE)
    except ModuleNotFoundError:
        return False
    draw_ranked_intervals_square(canvas, records)
    canvas.save(path)
    return True


def main() -> None:
    args = parse_args()
    raw_rows, records = load_stats(args.stats_csv)
    if not records:
        raise ValueError("No RNA quantile boxplot summaries were found.")

    output_dir = args.output_dir or args.stats_csv.parent / "nature_visualizations"
    output_dir.mkdir(parents=True, exist_ok=True)

    svg_path = output_dir / f"{args.basename}.svg"
    png_path = output_dir / f"{args.basename}.png"
    panel_a_svg_path = output_dir / f"{args.basename}_panel_a_square.svg"
    panel_a_png_path = output_dir / f"{args.basename}_panel_a_square.png"
    write_svg(svg_path, raw_rows, records)
    png_written = write_png(png_path, raw_rows, records)
    write_panel_a_square_svg(panel_a_svg_path, records)
    panel_a_png_written = write_panel_a_square_png(panel_a_png_path, records)

    iqr_values = [record.iqr for record in records]
    range_values = [record.full_range for record in records]
    print(f"Rows in CSV: {len(raw_rows):,}")
    print(f"Unique SNVs: {len({row.get('input_snv', '') for row in raw_rows if row.get('input_snv', '')}):,}")
    print(f"RNA summaries used: {len(records):,}")
    print(f"Median RNA IQR: {quantile(iqr_values, 0.50):.3f}")
    print(f"SNVs with RNA IQR > 0.30: {pct(sum(value > 0.30 for value in iqr_values) / len(iqr_values))}")
    print(f"Median RNA min-max range: {quantile(range_values, 0.50):.3f}")
    print(f"SNVs with RNA min-max range > 0.90: {pct(sum(value > 0.90 for value in range_values) / len(range_values))}")
    print(f"Wrote SVG: {svg_path}")
    if png_written:
        print(f"Wrote PNG: {png_path}")
    else:
        print("PNG preview skipped because PIL is not available.")
    print(f"Wrote square panel A SVG: {panel_a_svg_path}")
    if panel_a_png_written:
        print(f"Wrote square panel A PNG: {panel_a_png_path}")
    else:
        print("Square panel A PNG preview skipped because PIL is not available.")


if __name__ == "__main__":
    main()
