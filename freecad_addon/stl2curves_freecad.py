"""stl2curves for FreeCAD: turn mesh objects (or STL / 3MF files) into solids with true
curved surfaces.

The conversion runs as a separate Python process, not inside FreeCAD. FreeCAD's own
OpenCascade (its Part module) is not the build stl2curves is written for (cadquery-ocp),
and loading both into one process risks clashing libraries; a separate process also
keeps FreeCAD responsive and lets Stop end the conversion (and its worker processes)
outright. By default that process is FreeCAD's own Python, with stl2curves' dependencies
installed on first use into a folder of their own (FreeCAD's packages are left alone);
any other Python that has them can be chosen instead. The stl2curves package itself is
the one in this add-on's folder (the add-on is the whole repository), so the add-on and
the converter are always the same version.

A mesh is converted in its own coordinates (the frame its holes and faces are usually
lined up with, which stl2curves snaps axes to) and the solid given the mesh's placement.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback

import FreeCAD as App
import FreeCADGui as Gui
from PySide import QtCore, QtGui, QtWidgets

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                  # the repository: holds the stl2curves package
ICON = os.path.join(HERE, "Resources", "icons", "stl2curves.svg")
PARAMS = "User parameter:BaseApp/Preferences/Mod/stl2curves"
COMMAND = "Stl2Curves_Convert"
DEPENDENCIES = ["cadquery-ocp>=8.0", "numpy", "scipy"]
DISK_MB = 750                                 # what DEPENDENCIES take up once installed
NO_WINDOW = 0x08000000 if os.name == "nt" else 0   # CREATE_NO_WINDOW: no console flashing up


def _params():
    return App.ParamGet(PARAMS)


def private_packages():
    """Where stl2curves' dependencies go for FreeCAD's own Python (one folder per Python
    version: compiled packages only load in the version they were built for)."""
    return os.path.join(App.getUserAppDataDir(), "stl2curves",
                        f"py{sys.version_info[0]}{sys.version_info[1]}")


def freecad_python():
    """FreeCAD's bundled Python (sys.executable is FreeCAD itself inside the GUI)."""
    names = ("python.exe",) if os.name == "nt" else ("python3", "python")
    home = App.getHomePath()
    for folder in (os.path.dirname(sys.executable), os.path.join(home, "bin"), home):
        for name in names:
            path = os.path.join(folder, name)
            if os.path.isfile(path):
                return path
    return None


def interpreter(custom):
    """(python, PYTHONPATH entries, isolated): FreeCAD's Python with the private packages,
    or the chosen one with its own."""
    if custom:
        return shutil.which(custom) or custom, [ROOT], False
    return freecad_python(), [private_packages(), ROOT], True


def environment(paths, isolated):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(paths)
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    if isolated:
        env["PYTHONNOUSERSITE"] = "1"     # nothing from the user's own site-packages
    else:
        env.pop("PYTHONHOME", None)       # FreeCAD's, if set, would point another Python astray
    return env


def missing(python, paths, isolated):
    """Why that Python can't run stl2curves ('' if it can)."""
    try:
        r = subprocess.run([python, "-c", "import OCP, numpy, scipy, stl2curves"], env=environment(paths, isolated),
                           capture_output=True, text=True, timeout=180, creationflags=NO_WINDOW)
    except (OSError, subprocess.SubprocessError) as e:
        return str(e)
    if r.returncode == 0:
        return ""
    lines = [s for s in (r.stderr or r.stdout).splitlines() if s.strip()]
    return lines[-1] if lines else f"exit code {r.returncode}"


def _guarded(slot):
    """A slot that reports its own errors in the dialog (an exception escaping a Qt slot
    only reaches FreeCAD's console and leaves the dialog stuck half-way)."""
    takes = slot.__code__.co_argcount - 1     # Qt passes extras (clicked's "checked")
    def run(self, *args, **kwargs):
        try:
            return slot(self, *args[:takes], **kwargs)
        except Exception:
            self.say(traceback.format_exc())
            self._idle()
            self.stage = "idle"
    run.__name__ = slot.__name__
    return run


def _remove(folder, tries=5):
    """Delete a folder, retrying for a while: worker processes of a stopped conversion
    can hold their files open for a moment after it ends."""
    shutil.rmtree(folder, ignore_errors=True)
    if os.path.exists(folder) and tries > 1:
        QtCore.QTimer.singleShot(3000, lambda: _remove(folder, tries - 1))


def _safe(text):
    return re.sub(r"[^\w\-]+", "_", text).strip("_")[:40] or "mesh"


class ConvertDialog(QtWidgets.QDialog):
    """Options, then the conversion's output as it runs; the results go into the document."""

    def __init__(self, objects, files, doc_name):
        super().__init__(Gui.getMainWindow())
        self.setWindowTitle("Mesh to Curved Solid (stl2curves)")
        self.setWindowIcon(QtGui.QIcon(ICON))
        self.objects = [(o.Document.Name, o.Name) for o in objects]
        self.files = list(files)
        self.doc_name = doc_name
        self.proc = None
        self.stage = "idle"
        self.tmp = None
        self.created = []
        p = _params()

        if objects:
            what = ", ".join(o.Label for o in objects[:4]) + (f" and {len(objects) - 4} more" if len(objects) > 4 else "")
            source = f"Converts {'the mesh' if len(objects) == 1 else f'{len(objects)} meshes'}: {what}"
        else:
            source = "Converts " + ", ".join(os.path.basename(f) for f in files[:4]) + \
                     (f" and {len(files) - 4} more" if len(files) > 4 else "")
        top = QtWidgets.QLabel(source)
        top.setWordWrap(True)

        options = QtWidgets.QGroupBox("Options")
        form = QtWidgets.QFormLayout(options)
        self.blends = QtWidgets.QCheckBox("Smooth freeform faces where no simple surface fits")
        self.blends.setChecked(p.GetBool("Blends", True))
        self.true_size = QtWidgets.QCheckBox("Rebuild at the apparent design size (round radii)")
        self.true_size.setToolTip("Undo an overall scale such as a 99% slicer scale or an inch/cm export,\n"
                                  "with radii snapped to round values")
        self.true_size.setChecked(p.GetBool("TrueSize", False))
        self.hide = QtWidgets.QCheckBox("Hide the meshes afterwards")
        self.hide.setChecked(p.GetBool("HideMeshes", True))
        self.hide.setVisible(bool(objects))
        self.time_limit = QtWidgets.QSpinBox()
        self.time_limit.setRange(0, 24 * 3600)
        self.time_limit.setSingleStep(60)
        self.time_limit.setSuffix(" s")
        self.time_limit.setSpecialValueText("no limit")
        self.time_limit.setValue(p.GetInt("TimeLimit", 600))
        self.time_limit.setToolTip("After about this long, stop fitting further freeform faces and hunting\n"
                                   "down patches that spoil the solid, and keep what checks out")
        self.python = QtWidgets.QLineEdit(p.GetString("Python", ""))
        self.python.setPlaceholderText("FreeCAD's own Python")
        self.python.setToolTip("Leave empty to use FreeCAD's own Python, with stl2curves' dependencies\n"
                               f"installed on first use into {private_packages()}.\n"
                               "Or choose any Python 3.11+ that has cadquery-ocp, numpy and scipy.")
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self._browse_python)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.python)
        row.addWidget(browse)
        form.addRow(self.blends)
        form.addRow(self.true_size)
        form.addRow(self.hide)
        form.addRow("Time limit:", self.time_limit)
        form.addRow("Python:", row)

        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.log.setMinimumHeight(220)
        self.log.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.busy = QtWidgets.QProgressBar()
        self.busy.setRange(0, 0)
        self.busy.setVisible(False)
        self.clock = QtCore.QElapsedTimer()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._tick)

        self.convert_button = QtWidgets.QPushButton("Convert")
        self.convert_button.setDefault(True)
        self.convert_button.clicked.connect(self.start)
        self.stop_button = QtWidgets.QPushButton("Stop")
        self.stop_button.clicked.connect(self.stop)
        self.stop_button.setVisible(False)
        self.close_button = QtWidgets.QPushButton("Close")
        self.close_button.clicked.connect(self.close)
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.busy, 1)
        buttons.addWidget(self.convert_button)
        buttons.addWidget(self.stop_button)
        buttons.addWidget(self.close_button)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(top)
        layout.addWidget(options)
        layout.addWidget(self.log, 1)
        layout.addLayout(buttons)
        self.resize(720, 560)

    # ---- running ----------------------------------------------------------------------

    def running(self):
        return self.proc is not None and self.proc.state() != QtCore.QProcess.NotRunning

    def say(self, text):
        self.log.appendPlainText(text.rstrip("\n"))
        self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())

    @_guarded
    def start(self):
        p = _params()
        p.SetBool("Blends", self.blends.isChecked())
        p.SetBool("TrueSize", self.true_size.isChecked())
        p.SetBool("HideMeshes", self.hide.isChecked())
        p.SetInt("TimeLimit", self.time_limit.value())
        p.SetString("Python", self.python.text().strip())
        python, paths, isolated = interpreter(self.python.text().strip())
        if not python:
            self.say("Can't find FreeCAD's Python: choose a Python under Options.")
            return
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            problem = missing(python, paths, isolated)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        if problem and isolated:
            answer = QtWidgets.QMessageBox.question(
                self, "Install stl2curves' dependencies?",
                "stl2curves needs OpenCascade (cadquery-ocp), numpy and scipy, which FreeCAD's own "
                f"Python doesn't have. Download them now (about {DISK_MB} MB on disk)?\n\n"
                f"They go into {private_packages()}; FreeCAD's own packages are not changed.")
            if answer == QtWidgets.QMessageBox.Yes:
                self._install(python)
            return
        if problem:
            self.say(f"{python} can't run stl2curves: {problem}\n"
                     "Install its dependencies there (pip install cadquery-ocp numpy scipy), or clear "
                     "the Python field to use FreeCAD's own.")
            return
        self._convert(python, paths, isolated)

    def _install(self, python):
        self.say(f"Installing {', '.join(DEPENDENCIES)} into {private_packages()} …")
        os.makedirs(private_packages(), exist_ok=True)
        args = ["-m", "pip", "install", "--disable-pip-version-check", "--upgrade",
                "--target", private_packages()] + DEPENDENCIES
        env = environment([], True)
        env.pop("PYTHONPATH")
        self._run("install", python, args, env)

    def _convert(self, python, paths, isolated):
        self.tmp = tempfile.mkdtemp(prefix="stl2curves_")
        out = os.path.join(self.tmp, "out")
        self.jobs = []                        # (document, mesh name, label, placement, STEP file)
        inputs = []
        for i, (doc_name, name) in enumerate(self.objects):
            doc = App.listDocuments().get(doc_name)
            o = doc and doc.getObject(name)
            if o is None:
                continue
            stem = f"{_safe(o.Label)}_{i + 1}"
            mesh = o.Mesh.copy()
            mesh.Placement = App.Placement()  # its own coordinates; the solid gets its placement
            path = os.path.join(self.tmp, stem + ".stl")
            mesh.write(path)
            inputs.append(path)
            self.jobs.append((doc_name, name, o.Label, o.getGlobalPlacement(), os.path.join(out, stem + ".step")))
        inputs += self.files
        if not inputs:
            self.say("Nothing left to convert (the meshes were deleted).")
            return
        args = ["-m", "stl2curves"] + inputs + ["--out", out, "--time-limit", str(self.time_limit.value())]
        if not self.blends.isChecked():
            args.append("--no-blends")
        if self.true_size.isChecked():
            args.append("--true-size")
        self.out = out
        env = environment(paths, isolated)
        # its temporary files (the mesh it shares with its workers, scratch folders) go
        # under ours: a stopped conversion can't tidy up after itself
        temp = os.path.join(self.tmp, "temp")
        os.makedirs(temp)
        env["TMP"] = env["TEMP"] = env["TMPDIR"] = temp
        self._run("convert", python, args, env)

    def _run(self, stage, python, args, env):
        self.stage = stage
        self.stopped = False
        self.pending = b""
        self.proc = QtCore.QProcess(self)
        self.proc.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        qenv = QtCore.QProcessEnvironment()
        for k, v in env.items():
            qenv.insert(k, v)
        self.proc.setProcessEnvironment(qenv)
        self.proc.setWorkingDirectory(self.tmp or os.path.expanduser("~"))
        self.proc.readyReadStandardOutput.connect(self._output)
        self.proc.finished.connect(self._finished)
        self.proc.errorOccurred.connect(self._error)
        self.convert_button.setEnabled(False)
        self.stop_button.setVisible(True)
        self.busy.setVisible(True)
        self.clock.start()
        self.timer.start(1000)
        self.proc.start(python, args)

    @_guarded
    def _output(self, flush=False):
        # whole lines only: a chunk can end mid-line (or mid-character)
        self.pending += bytes(self.proc.readAllStandardOutput())
        cut = len(self.pending) if flush else self.pending.rfind(b"\n") + 1
        if cut <= 0:
            return
        text, self.pending = self.pending[:cut].decode("utf-8", "replace"), self.pending[cut:]
        text = text.replace("\r\n", "\n").replace("\r", "\n")
        if text.strip():
            self.say(text)

    def _tick(self):
        s = self.clock.elapsed() // 1000
        self.busy.setFormat(f"{s // 60}:{s % 60:02d}")
        self.busy.setTextVisible(True)

    def _idle(self):
        self.timer.stop()
        self.busy.setVisible(False)
        self.stop_button.setVisible(False)
        self.convert_button.setEnabled(True)

    @_guarded
    def _error(self, error):
        if error == QtCore.QProcess.FailedToStart:
            self.say(f"Couldn't start {self.proc.program()}: {self.proc.errorString()}")
            self._idle()
            self.stage = "idle"

    @_guarded
    def _finished(self, code, status):
        stage, self.stage = self.stage, "idle"
        self._output(flush=True)
        self._idle()
        crashed = status != QtCore.QProcess.NormalExit
        if stage == "install":
            if code or crashed:
                self.say("Installing the dependencies failed (see above).")
                return
            self.say("Dependencies installed.\n")
            self.start()
            return
        if stage == "convert":
            if crashed:
                self.say("Stopped." if self.stopped else "The conversion crashed.")
            elif code:
                self.say(f"The conversion ended with an error (exit code {code}).")
            self._take_results()
            _remove(self.tmp)
            self.tmp = None

    def _take_results(self):
        claimed = set()
        new = []
        docs = set()
        for doc_name, name, label, placement, step in self.jobs:
            claimed.add(os.path.normcase(step))
            doc = App.listDocuments().get(doc_name)
            if doc is None or not os.path.isfile(step):
                continue
            new.append(self._add(doc, step, f"{label} solid", placement))
            docs.add(doc)
            mesh = doc.getObject(name)
            if mesh is not None and self.hide.isChecked() and mesh.ViewObject is not None:
                mesh.ViewObject.Visibility = False
        others = sorted(f for f in (os.listdir(self.out) if os.path.isdir(self.out) else [])
                        if f.lower().endswith(".step") and os.path.normcase(os.path.join(self.out, f)) not in claimed)
        if others:
            doc = App.listDocuments().get(self.doc_name) if self.doc_name else None
            doc = doc or App.newDocument()
            for f in others:
                new.append(self._add(doc, os.path.join(self.out, f), os.path.splitext(f)[0], App.Placement()))
                docs.add(doc)
        for doc in docs:
            doc.recompute()
        if new:
            Gui.Selection.clearSelection()
            for o in new:
                Gui.Selection.addSelection(o)
            self.say("\nAdded " + ", ".join(o.Label for o in new) + ".")
        self.created = new

    def _add(self, doc, step, label, placement):
        import Part
        shape = Part.Shape()
        shape.read(step)
        doc.openTransaction("stl2curves")
        o = doc.addObject("Part::Feature", "CurvedSolid")
        o.Shape = shape
        o.Label = label
        o.Placement = placement
        doc.commitTransaction()
        return o

    def stop(self):
        if self.running():
            self.stopped = True
            self.say("Stopping …")
            self.proc.kill()      # its worker processes go with it (stl2curves.workers)

    def _browse_python(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Python with stl2curves' dependencies",
                                                        self.python.text() or os.path.expanduser("~"),
                                                        "Python (python*.exe python*);;All files (*)")
        if path:
            self.python.setText(path)

    def closeEvent(self, event):
        if self.running():
            if QtWidgets.QMessageBox.question(self, "stl2curves", "Stop the conversion?") != QtWidgets.QMessageBox.Yes:
                event.ignore()
                return
            self.stop()
            self.proc.waitForFinished(5000)
        _open.discard(self)
        event.accept()


_open = set()       # dialogs shown (modeless: kept here so Python doesn't collect them)


def open_dialog(objects=None, files=None):
    if objects is None and files is None:
        objects = [o for o in Gui.Selection.getSelection() if o.isDerivedFrom("Mesh::Feature")]
        files = []
        if not objects:
            p = _params()
            files, _ = QtWidgets.QFileDialog.getOpenFileNames(
                Gui.getMainWindow(), "Meshes to convert (or select mesh objects first)",
                p.GetString("LastFolder", os.path.expanduser("~")),
                "Meshes (*.stl *.STL *.3mf *.3MF);;All files (*)")
            if not files:
                return None
            p.SetString("LastFolder", os.path.dirname(files[0]))
    objects, files = objects or [], files or []
    doc = objects[0].Document.Name if objects else (App.ActiveDocument.Name if App.ActiveDocument else None)
    dialog = ConvertDialog(objects, files, doc)
    _open.add(dialog)
    dialog.show()
    return dialog


class ConvertCommand:
    def GetResources(self):
        return {"Pixmap": ICON,
                "MenuText": "Mesh to Curved Solid (stl2curves)…",
                "ToolTip": "Turn the selected meshes (or STL / 3MF files, if none is selected) into solids\n"
                           "with true curved surfaces: holes, fillets, chamfers, threads and smooth areas\n"
                           "become real cylinders, cones, tori, helices and B-spline faces"}

    def IsActive(self):
        return True

    def Activated(self):
        open_dialog()


class MenuManipulator:
    """Puts the command next to FreeCAD's own mesh-to-shape commands in every workbench
    that has them (Part, Mesh), and on their toolbars."""

    def modifyMenuBar(self):
        return [{"insert": COMMAND, "menuItem": "Part_ShapeFromMesh", "after": ""},
                {"insert": COMMAND, "menuItem": "Mesh_FromPartShape", "after": ""}]

    def modifyToolBars(self):
        return [{"append": COMMAND, "toolBar": "Part Tools"},
                {"append": COMMAND, "toolBar": "Mesh Tools"}]


def setup():
    Gui.addCommand(COMMAND, ConvertCommand())
    Gui.addWorkbenchManipulator(MenuManipulator())
