module Library

  # REQUIRED METHODS (CALLED BY DRIVER)

  def lib_build(env)
    # Construct the command to execute in a subshell to perform the build.
    cmd="cd #{env.build.dir} && make"
    # Execute external command via ext() (defined in ts.rb).
    ext(cmd,{:msg=>"Build failed, see #{logfile}"})
  end

  def lib_build_post(env,output)
    # Move the executables into a directory of their own.
    bindir=File.join(env.build.dir,"bin")
    FileUtils.mkdir_p(bindir)
    [
      "Do_Smv2_Solver.exe",
      "GearProc/GearProc.exe"
    ].each do |x|
      FileUtils.mv(File.join(env.build.dir,x),bindir)
      logd "Moved #{x} -> #{bindir}"
    end
    # Return the name of the bin dir to be copied for each run that requires it.
    bindir
  end

  def lib_build_prep(env)
    # Construct the name of the build directory and store it in the env.build
    # structure for later reference. The value of env.build.root is supplied
    # internally by the test suite; the value of env.run.build is supplied by
    # the run config.
    env.build.dir=File.join(env.build.root,env.run.build)
    # Construct the path to the source files. Wrapping in valid_dir() (defined
    # in ts.rb) ensures that it actually already exists.
    srcdir=valid_dir(File.join("..","src"))
    # Copy the source files, recursively, into the build directory.
    FileUtils.cp_r(srcdir,env.build.dir)
    logd "Copied #{srcdir} -> #{env.build.dir}"
    # Link the correct top-level makefile.
    src=env.build.makefile_top
    dst=File.join(env.build.dir,"GNUmakefile")
    FileUtils.ln_sf(src,dst)
    logd "Linked top-level makefile: #{dst} -> #{src}"
    # Link the correct GearProc makefile.
    src=env.build.makefile_gearproc
    dst=File.join(env.build.dir,"GearProc","GNUmakefile")
    FileUtils.ln_sf(src,dst)
    logd "Linked GearProc makefile: #{dst} -> #{src}"
  end

  def lib_outfiles(env,path)
    [
      [path,"physproc_exit.proc0017"],
      [path,"smv2chem1_exit.proc0017"],
      [path,"smv2chem2_exit.proc0017"]
    ]
  end

  def lib_prep_data(env)
    logd "No data-prep needed."
  end

  def lib_prep_job(env,rundir)
    # Copy executable dir into run directory. The value of env.build.runfiles
    # is provided internally by the test suite.
    FileUtils.cp_r(env.build.runfiles,rundir)
    logd "Copied #{env.build.runfiles} -> #{rundir}"
    # Since the executables are in a subdirectory of the run directory, update
    # rundir to reflect this.
    rundir=File.join(rundir,File.basename(env.build.runfiles))
    # Link data.
    datadir=valid_dir(File.join("..","data"))
    [
      "physproc_entry.proc0017",
      "smv2chem1_entry.proc0017",
      "smv2chem2_entry.proc0017"
    ].each { |x| FileUtils.ln_sf(File.join(datadir,x),rundir) }
    # Return rundir (i.e. where to perform the run).
    rundir
  end

  def lib_queue_del_cmd(env)
    nil
  end

  def lib_re_str_success(env)
    "Exiting doSmv2Solver"
  end

  def lib_run_job(env,rundir)
    # Construct the command to execute in a subshell to perform the run.
    cmd="cd #{rundir} && ./Do_Smv2_Solver.exe > stdout"
    # Execute external command via ext() (defined in ts.rb).
    ext(cmd,{:msg=>"Run failed, see #{logfile}"})
    # Return the path to the run's 'stdout' file.
    File.join(rundir,"stdout")
  end

  # CUSTOM METHODS (NOT CALLED BY DRIVER)

end
