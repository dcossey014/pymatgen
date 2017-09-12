# coding: utf-8

from __future__ import unicode_literals

from fireworks import LaunchPad
from custodian.ansible.actions import FileActions, DictActions
from custodian.ansible.interpreter import Modder
from pymatgen.io.bgw import BgwInput


class BgwModder(Modder):
    def __init__(self, actions=None, strict=True, bgw=None):
        """
        Initializes a Modder for VaspInput sets

        Args:
            actions ([Action]): A sequence of supported actions. See
                :mod:`custodian.ansible.actions`. Default is None,
                which means DictActions and FileActions are supported.
            strict (bool): Indicating whether to use strict mode. In non-strict
                mode, unsupported actions are simply ignored without any
                errors raised. In strict mode, if an unsupported action is
                supplied, a ValueError is raised. Defaults to True.
            bgw (BgwInput): A BgwInput object from the current directory.
                Initialized automatically if not passed (but passing it will
                avoid having to reparse the directory).
        """
        self.bgw = bgw or BgwInput.from_directory('.')
        actions = actions or [FileActions, DictActions]
        super(BgwModder, self).__init__(actions, strict)

    def apply_actions(self, actions):
        """
        Applies a list of actions to the Vasp Input Set and rewrites modified
        files.
        Args:
            actions [dict]: A list of actions of the form {'file': filename,
                'action': moddermodification} or {'dict': bgw_input_obj,
                'action': moddermodification}
        """
        for a in actions:
            if "launchpad" in a:
                lp = a['launchpad']
                fw = lp.get_fw_by_id(a['fw_id'])

                # Set Firework state to Ready so we can modify it
                if fw.state != "READY":
                    lp.rerun_fw(a['fw_id'])

                lp.update_spec([a['fw_id']], a['action'])
            elif "dict" in a:
                k = a["dict"]
                modified.append(k)
                self.bgw = self.modify_object(a["action"], self.bgw)
                self.bgw.to_file(self.bgw.filename)
            elif "file" in a:
                self.modify(a["action"], a["file"])
            else:
                raise ValueError("Unrecognized format: {}".format(a))
